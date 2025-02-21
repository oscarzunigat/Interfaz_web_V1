from flask import Blueprint, render_template, request

pagina2_bp = Blueprint('pagina2', __name__)

@pagina2_bp.route('/pagina2')
def pagina2():
    return render_template('pagina2.html')

@pagina2_bp.route('/suma2', methods=['POST'])
def suma2():
    num1 = request.form.get('num1', type=int)
    num2 = request.form.get('num2', type=int)
    resultado = num1 + num2
    return f'La suma de {num1} y {num2} es {resultado}.'
