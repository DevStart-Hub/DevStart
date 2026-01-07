import sys
import random

print("🧪 Testing PsychoPy installation...")
print("This will open a window with a countdown and celebration message!")

# --- STEP 1: CRITICAL IMPORTS (Visual & Core) ---
try:
    # We import event separately as it's often needed for keypresses
    from psychopy import visual, core, event
    print("✅ Critical modules (visual, core, event) imported successfully.")
except ImportError as e:
    print(f"\n❌ CRITICAL ERROR: Could not import basic PsychoPy modules.")
    print(f"Error details: {e}")
    sys.exit(1)

# --- STEP 2: OPTIONAL SOUND IMPORT ---
# We try to import and initialize sound separately. 
# If this fails, we just disable sound but keep running.
sound_working = False
beep_sound = None
success_sound = None

try:
    from psychopy import sound
    print("✅ Sound module imported.")
    
    # Try to initialize a dummy sound to test the audio backend
    # This catches errors where the module imports but the driver fails
    test_sound = sound.Sound('C', secs=0.1)
    
    # If we got here, sound is safe to use
    sound_working = True
    print("✅ Audio backend initialized successfully.")
    
except Exception as e:
    print(f"⚠️  SOUND WARNING: Sound system encountered an error.")
    print(f"   Error details: {e}")
    print("   -> Continuing without sound (Visuals will still work!)")
    sound_working = False

# --- STEP 3: THE EXPERIMENT LOOP ---
try:
    # Create a window
    win = visual.Window(
        size=(800, 600),
        color='lightgray',
        units='pix',
        fullscr=False,
        allowGUI=True
    )
    
    # Create title text
    title_text = visual.TextStim(
        win, text='PsychoPy Test', color='darkblue',
        height=60, pos=(0, 200), bold=True
    )
    
    # Create the countdown text stimulus
    countdown_text = visual.TextStim(
        win, text='', color='red',
        height=120, pos=(0, 0), bold=True
    )
    
    # Create instruction text
    instruction_text = visual.TextStim(
        win, text='Get ready...', color='black',
        height=30, pos=(0, -150)
    )
    
    # Initialize Sounds (Only if the import/test above worked)
    if sound_working:
        try:
            beep_sound = sound.Sound('C', secs=0.2, hamming=True)
            success_sound = sound.Sound('A', secs=0.5, hamming=True)
        except:
            # Fallback if sound fails during creation
            sound_working = False
    
    # Show title and instructions
    title_text.draw()
    instruction_text.draw()
    win.flip()
    core.wait(2)
    
    # Show countdown from 3 to 1
    for i in range(3, 0, -1):
        title_text.draw()
        
        countdown_text.text = str(i)
        
        if i == 3: countdown_text.color = 'red'
        elif i == 2: countdown_text.color = 'orange'
        else: countdown_text.color = 'green'
        
        countdown_text.draw()
        win.flip()
        
        if sound_working:
            beep_sound.play()
        
        core.wait(1)
    
    # Celebration Confetti
    circles = []
    colors = ['red', 'blue', 'green', 'yellow', 'purple', 'orange', 'pink']
    
    for i in range(20):
        circle = visual.Circle(
            win,
            radius=random.randint(10, 30),
            pos=(random.randint(-300, 300), random.randint(-200, 200)),
            fillColor=random.choice(colors),
            lineColor=None
        )
        circles.append(circle)
    
    # Success Message
    msg = '🎉 FANTASTIC! 🎉\n\nVisuals are working perfectly!'
    if not sound_working:
        msg += '\n(Sound is disabled, but that is OK)'
    
    success_message = visual.TextStim(
        win, text=msg, color='darkgreen',
        height=30, pos=(0, 0), bold=True
    )
    
    # Info Text
    status_msg = '✅ Visual stimuli: Working\n✅ Timing: Working\n'
    if sound_working:
        status_msg += '✅ Sound: Working'
    else:
        status_msg += '⚠️ Sound: Failed (Software/Driver issue)'
        
    info_text = visual.TextStim(
        win, text=status_msg, color='black',
        height=25, pos=(0, -200)
    )
    
    # Animation Loop
    for frame in range(180):
        for circle in circles:
            circle.draw()
        
        success_message.draw()
        info_text.draw()
        win.flip()
        
        if frame == 0 and sound_working:
            success_sound.play()
            
        core.wait(0.0167)
    
    # Final cleanup instruction
    final_text = visual.TextStim(
        win, text='Press any key to close...', color='gray',
        height=25, pos=(0, -250)
    )
    final_text.draw()
    success_message.draw()
    info_text.draw()
    win.flip()
    
    # Wait for key press
    event.waitKeys()
    
    # --- FINAL REPORT ---
    print("\n" + "="*40)
    print("       TEST RESULTS       ")
    print("="*40)
    print("✅ Visual System: WORKING")
    print("✅ Event/Timing:  WORKING")
    
    if sound_working:
        print("✅ Sound System:  WORKING")
        print("\n🎉 COMPLETELY SUCCESSFUL! You are 100% ready.")
    else:
        print("⚠️  Sound System:  NOT WORKING")
        print("\n👉 DIAGNOSIS: The sound library failed to load.")
        print("   HOWEVER: This is often just a driver issue.")
        print("   FOR THE WORKSHOP: Don't worry! You can likely proceed.")
        print("   (Most beginner tutorials focus on visual stimuli anyway).")
    print("="*40 + "\n")

except Exception as e:
    print(f"\n❌ ERROR: Something went wrong during the test.")
    print(f"Error message: {e}")

finally:
    try:
        win.close()
    except:
        pass
    core.quit()