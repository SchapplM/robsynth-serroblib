% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:18:43
% EndTime: 2018-11-23 16:18:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->42)
t27 = sin(qJ(4));
t46 = t27 * pkin(4);
t29 = cos(qJ(4));
t11 = t29 * pkin(4) + pkin(3);
t20 = pkin(10) + qJ(3);
t14 = cos(t20);
t28 = sin(qJ(1));
t45 = t28 * t14;
t21 = qJ(4) + pkin(11);
t15 = cos(t21);
t44 = t28 * t15;
t43 = t28 * t27;
t42 = t28 * t29;
t30 = cos(qJ(1));
t41 = t30 * t14;
t40 = t30 * t15;
t39 = t30 * t27;
t38 = t30 * t29;
t25 = -qJ(5) - pkin(8);
t22 = pkin(6) + 0;
t24 = cos(pkin(10));
t9 = t24 * pkin(2) + pkin(1);
t37 = t30 * t9 + 0;
t23 = sin(pkin(10));
t36 = t23 * pkin(2) + t22;
t26 = -pkin(7) - qJ(2);
t35 = t30 * t26 + t28 * t9 + 0;
t12 = sin(t20);
t34 = pkin(3) * t14 + pkin(8) * t12;
t1 = pkin(5) * t15 + t11;
t19 = -pkin(9) + t25;
t33 = t1 * t14 - t12 * t19;
t32 = t11 * t14 - t12 * t25;
t31 = -t28 * t26 + t37;
t16 = qJ(6) + t21;
t13 = sin(t21);
t8 = cos(t16);
t7 = sin(t16);
t6 = t30 * t12;
t5 = t28 * t12;
t2 = pkin(5) * t13 + t46;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t30 * t24, -t30 * t23, t28, t30 * pkin(1) + t28 * qJ(2) + 0; t28 * t24, -t28 * t23, -t30, t28 * pkin(1) - t30 * qJ(2) + 0; t23, t24, 0, t22; 0, 0, 0, 1; t41, -t6, t28, t31; t45, -t5, -t30, t35; t12, t14, 0, t36; 0, 0, 0, 1; t14 * t38 + t43, -t14 * t39 + t42, t6, t30 * t34 + t31; t14 * t42 - t39, -t14 * t43 - t38, t5, t28 * t34 + t35; t12 * t29, -t12 * t27, -t14, t12 * pkin(3) - t14 * pkin(8) + t36; 0, 0, 0, 1; t28 * t13 + t14 * t40, -t13 * t41 + t44, t6, t32 * t30 + (-t26 + t46) * t28 + t37; -t30 * t13 + t14 * t44, -t13 * t45 - t40, t5, -pkin(4) * t39 + t28 * t32 + t35; t12 * t15, -t12 * t13, -t14, t12 * t11 + t14 * t25 + t36; 0, 0, 0, 1; t28 * t7 + t8 * t41, t28 * t8 - t7 * t41, t6, t33 * t30 + (t2 - t26) * t28 + t37; -t30 * t7 + t8 * t45, -t30 * t8 - t7 * t45, t5, -t30 * t2 + t28 * t33 + t35; t12 * t8, -t12 * t7, -t14, t12 * t1 + t14 * t19 + t36; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
