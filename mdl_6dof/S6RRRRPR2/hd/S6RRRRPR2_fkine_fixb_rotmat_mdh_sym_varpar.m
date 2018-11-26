% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:12:20
% EndTime: 2018-11-23 18:12:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->59), mult. (112->64), div. (0->0), fcn. (166->12), ass. (0->40)
t30 = -pkin(8) - pkin(7);
t28 = cos(qJ(2));
t11 = t28 * pkin(2) + pkin(1);
t19 = pkin(11) + qJ(6);
t12 = sin(t19);
t27 = sin(qJ(1));
t45 = t27 * t12;
t13 = cos(t19);
t44 = t27 * t13;
t23 = sin(pkin(11));
t43 = t27 * t23;
t24 = cos(pkin(11));
t42 = t27 * t24;
t29 = cos(qJ(1));
t41 = t29 * t12;
t40 = t29 * t13;
t39 = t29 * t23;
t38 = t29 * t24;
t22 = qJ(2) + qJ(3);
t20 = pkin(6) + 0;
t15 = cos(t22);
t3 = pkin(3) * t15 + t11;
t37 = t29 * t3 + 0;
t21 = -pkin(9) + t30;
t36 = t29 * t21 + t27 * t3 + 0;
t26 = sin(qJ(2));
t35 = t26 * pkin(2) + t20;
t14 = sin(t22);
t34 = pkin(3) * t14 + t35;
t16 = qJ(4) + t22;
t10 = cos(t16);
t25 = -pkin(10) - qJ(5);
t7 = t24 * pkin(5) + pkin(4);
t9 = sin(t16);
t33 = t10 * t7 - t25 * t9;
t32 = pkin(4) * t10 + qJ(5) * t9;
t31 = -t27 * t21 + t37;
t5 = t29 * t9;
t4 = t27 * t9;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t29 * t28, -t29 * t26, t27, t29 * pkin(1) + t27 * pkin(7) + 0; t27 * t28, -t27 * t26, -t29, t27 * pkin(1) - t29 * pkin(7) + 0; t26, t28, 0, t20; 0, 0, 0, 1; t29 * t15, -t29 * t14, t27, t29 * t11 - t27 * t30 + 0; t27 * t15, -t27 * t14, -t29, t27 * t11 + t29 * t30 + 0; t14, t15, 0, t35; 0, 0, 0, 1; t29 * t10, -t5, t27, t31; t27 * t10, -t4, -t29, t36; t9, t10, 0, t34; 0, 0, 0, 1; t10 * t38 + t43, -t10 * t39 + t42, t5, t29 * t32 + t31; t10 * t42 - t39, -t10 * t43 - t38, t4, t27 * t32 + t36; t9 * t24, -t9 * t23, -t10, t9 * pkin(4) - t10 * qJ(5) + t34; 0, 0, 0, 1; t10 * t40 + t45, -t10 * t41 + t44, t5, t33 * t29 + (pkin(5) * t23 - t21) * t27 + t37; t10 * t44 - t41, -t10 * t45 - t40, t4, -pkin(5) * t39 + t27 * t33 + t36; t9 * t13, -t9 * t12, -t10, t10 * t25 + t9 * t7 + t34; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
