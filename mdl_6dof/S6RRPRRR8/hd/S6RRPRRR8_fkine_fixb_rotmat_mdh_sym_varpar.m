% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:25:34
% EndTime: 2018-11-23 17:25:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (206->73), mult. (164->84), div. (0->0), fcn. (228->12), ass. (0->38)
t27 = sin(pkin(11));
t19 = t27 * pkin(3);
t28 = cos(pkin(11));
t11 = t28 * pkin(3) + pkin(2);
t25 = pkin(11) + qJ(4);
t15 = sin(t25);
t4 = pkin(4) * t15 + t19;
t31 = sin(qJ(1));
t43 = t31 * t27;
t32 = cos(qJ(2));
t42 = t31 * t32;
t33 = cos(qJ(1));
t41 = t33 * t32;
t29 = -pkin(8) - qJ(3);
t26 = pkin(6) + 0;
t16 = cos(t25);
t3 = pkin(4) * t16 + t11;
t40 = t31 * pkin(1) + 0;
t24 = -pkin(9) + t29;
t39 = t33 * pkin(1) + t31 * pkin(7) + 0;
t17 = qJ(5) + t25;
t10 = cos(t17);
t1 = pkin(5) * t10 + t3;
t18 = -pkin(10) + t24;
t30 = sin(qJ(2));
t38 = t1 * t32 - t18 * t30;
t37 = -t24 * t30 + t3 * t32;
t36 = pkin(2) * t32 + qJ(3) * t30;
t35 = t11 * t32 - t29 * t30;
t34 = -t33 * pkin(7) + t40;
t14 = t33 * t30;
t13 = t31 * t30;
t12 = qJ(6) + t17;
t9 = sin(t17);
t6 = cos(t12);
t5 = sin(t12);
t2 = pkin(5) * t9 + t4;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t31, 0, 0; t31, t33, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t41, -t14, t31, t39; t42, -t13, -t33, t34; t30, t32, 0, t26; 0, 0, 0, 1; t28 * t41 + t43, -t27 * t41 + t31 * t28, t14, t36 * t33 + t39; -t33 * t27 + t28 * t42, -t27 * t42 - t33 * t28, t13, t36 * t31 + t34; t30 * t28, -t30 * t27, -t32, t30 * pkin(2) - t32 * qJ(3) + t26; 0, 0, 0, 1; t31 * t15 + t16 * t41, -t15 * t41 + t31 * t16, t14, pkin(3) * t43 + t35 * t33 + t39; -t33 * t15 + t16 * t42, -t15 * t42 - t33 * t16, t13 (-pkin(7) - t19) * t33 + t35 * t31 + t40; t30 * t16, -t30 * t15, -t32, t30 * t11 + t32 * t29 + t26; 0, 0, 0, 1; t10 * t41 + t31 * t9, t31 * t10 - t9 * t41, t14, t31 * t4 + t37 * t33 + t39; t10 * t42 - t33 * t9, -t33 * t10 - t9 * t42, t13 (-pkin(7) - t4) * t33 + t37 * t31 + t40; t30 * t10, -t30 * t9, -t32, t32 * t24 + t30 * t3 + t26; 0, 0, 0, 1; t31 * t5 + t6 * t41, t31 * t6 - t5 * t41, t14, t31 * t2 + t38 * t33 + t39; -t33 * t5 + t6 * t42, -t33 * t6 - t5 * t42, t13 (-pkin(7) - t2) * t33 + t38 * t31 + t40; t30 * t6, -t30 * t5, -t32, t30 * t1 + t32 * t18 + t26; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
