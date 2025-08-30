from psychopy import visual, core, sound
import sys

print("🧪 Testing PsychoPy installation...")
print("This will open a window with a countdown and celebration message!")

try:
    # Create a window with better settings
    win = visual.Window(
        size=(800, 600),
        color='lightgray',
        units='pix',
        fullscr=False,
        allowGUI=True
    )
    
    # Create title text
    title_text = visual.TextStim(
        win,
        text='PsychoPy Test',
        color='darkblue',
        height=60,
        pos=(0, 200),
        bold=True
    )
    
    # Create the countdown text stimulus
    countdown_text = visual.TextStim(
        win,
        text='',
        color='red',
        height=120,
        pos=(0, 0),
        bold=True
    )
    
    # Create instruction text
    instruction_text = visual.TextStim(
        win,
        text='Get ready...',
        color='black',
        height=30,
        pos=(0, -150)
    )
    
    # Create sounds (using built-in sounds that should work on most systems)
    try:
        # Countdown beep sound
        beep_sound = sound.Sound('C', secs=0.2, hamming=True)
        # Success sound (higher pitch)
        success_sound = sound.Sound('A', secs=0.5, hamming=True)
        sounds_available = True
        print("✅ Sound system working!")
    except Exception as sound_error:
        print(f"⚠️  Sound not available: {sound_error}")
        print("   (This is okay - continuing without sound)")
        sounds_available = False
    
    # Show title and instructions first
    title_text.draw()
    instruction_text.draw()
    win.flip()
    core.wait(2)
    
    # Show countdown from 3 to 1
    for i in range(3, 0, -1):
        # Clear and draw title
        title_text.draw()
        
        # Draw countdown number
        countdown_text.text = str(i)
        countdown_text.draw()
        
        # Add some visual flair - change color for each number
        if i == 3:
            countdown_text.color = 'red'
        elif i == 2:
            countdown_text.color = 'orange'
        else:
            countdown_text.color = 'green'
        
        win.flip()
        
        # Play countdown sound
        if sounds_available:
            beep_sound.play()
        
        core.wait(1)
    
    # Create celebration visual effects
    # Create multiple colored circles for confetti effect
    circles = []
    import random
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
    
    # Create the final message text stimulus
    success_message = visual.TextStim(
        win,
        text='🎉 FANTASTIC! 🎉\n\nYour PsychoPy installation\nis working perfectly!\n\nYou are ready for\nthe DevStart tutorials!',
        color='darkgreen',
        height=35,
        pos=(0, 0),
        bold=True
    )
    
    # Create additional info text
    info_text = visual.TextStim(
        win,
        text='✅ Visual stimuli: Working\n✅ Timing: Working\n' + 
             ('✅ Sound: Working' if sounds_available else '⚠️  Sound: Not available'),
        color='black',
        height=25,
        pos=(0, -200)
    )
    
    # Show celebration screen with effects
    for frame in range(180):  # 3 seconds at 60fps
        # Draw confetti circles
        for circle in circles:
            circle.draw()
        
        # Draw success message
        success_message.draw()
        info_text.draw()
        
        win.flip()
        
        # Play success sound once at the beginning
        if frame == 0 and sounds_available:
            success_sound.play()
        
        # Small delay to control frame rate
        core.wait(0.0167)  # ~60fps
    
    # Final instruction
    final_text = visual.TextStim(
        win,
        text='Press any key or wait 3 seconds to close...',
        color='gray',
        height=25,
        pos=(0, -250)
    )
    
    final_text.draw()
    success_message.draw()
    info_text.draw()
    win.flip()
    
    # Wait for key press or timeout
    keys = core.event.waitKeys(maxWait=3.0)
    
    print("\n🎉 SUCCESS! PsychoPy is working perfectly!")
    print("✅ Window creation: Working")
    print("✅ Text rendering: Working") 
    print("✅ Visual stimuli: Working")
    print("✅ Timing functions: Working")
    sound_status = 'Working' if sounds_available else 'Not available (but that\'s okay)'
    print(f"✅ Sound system: {sound_status}")
    print("\n🚀 You're ready for PsychoPy experiments!")
    
except Exception as e:
    print(f"\n❌ ERROR: PsychoPy test failed!")
    print(f"Error message: {e}")
    print("\n🔧 Troubleshooting tips:")
    print("1. Make sure PsychoPy is properly installed")
    print("2. Try: pip install psychopy")
    print("3. Check if your graphics drivers are up to date")
    print("4. Make sure you're not running this in a headless environment")
    print("5. Try restarting your Python environment")

finally:
    # Clean up
    try:
        win.close()
    except:
        pass
    core.quit()

print("\n--- PsychoPy test completed ---")