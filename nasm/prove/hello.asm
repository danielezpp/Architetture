section .data
    msg db "Hello, World!", 0xA  ; Stringa da stampare con newline
    len equ $ - msg              ; Calcolo della lunghezza della stringa

section .text
    global main                ; Entry point per il linker

main:
    ; syscall write (stdout)
    mov eax, 4                   ; syscall: write
    mov ebx, 1                   ; file descriptor: stdout
    mov ecx, msg                 ; Puntatore alla stringa
    mov edx, len                 ; Lunghezza della stringa
    int 0x80                     ; Interrompi per la syscall

    ; syscall exit
    mov eax, 1                   ; syscall: exit
    xor ebx, ebx                 ; Exit code 0
    int 0x80                     ; Interrompi per la syscall