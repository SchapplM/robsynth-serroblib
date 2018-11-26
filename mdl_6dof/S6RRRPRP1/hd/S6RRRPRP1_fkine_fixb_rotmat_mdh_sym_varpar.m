% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:41:04
% EndTime: 2018-11-23 17:41:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (190->54), mult. (112->54), div. (0->0), fcn. (166->10), ass. (0->39)
t32 = -pkin(8) - pkin(7);
t30 = cos(qJ(2));
t16 = t30 * pkin(2) + pkin(1);
t24 = qJ(2) + qJ(3);
t17 = pkin(10) + t24;
t12 = sin(t17);
t26 = sin(qJ(5));
t44 = t12 * t26;
t28 = sin(qJ(1));
t43 = t28 * t26;
t29 = cos(qJ(5));
t42 = t28 * t29;
t31 = cos(qJ(1));
t41 = t31 * t26;
t40 = t31 * t29;
t23 = pkin(6) + 0;
t19 = cos(t24);
t7 = pkin(3) * t19 + t16;
t39 = t31 * t7 + 0;
t27 = sin(qJ(2));
t38 = t27 * pkin(2) + t23;
t22 = -qJ(4) + t32;
t37 = t31 * t22 + t28 * t7 + 0;
t18 = sin(t24);
t36 = pkin(3) * t18 + t38;
t13 = cos(t17);
t35 = pkin(4) * t13 + pkin(9) * t12;
t15 = t29 * pkin(5) + pkin(4);
t25 = -qJ(6) - pkin(9);
t34 = -t12 * t25 + t13 * t15;
t33 = -t28 * t22 + t39;
t10 = t31 * t12;
t9 = t28 * t12;
t8 = t12 * t29;
t4 = t13 * t40 + t43;
t3 = -t13 * t41 + t42;
t2 = t13 * t42 - t41;
t1 = -t13 * t43 - t40;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t28, 0, 0; t28, t31, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t31 * t30, -t31 * t27, t28, t31 * pkin(1) + t28 * pkin(7) + 0; t28 * t30, -t28 * t27, -t31, t28 * pkin(1) - t31 * pkin(7) + 0; t27, t30, 0, t23; 0, 0, 0, 1; t31 * t19, -t31 * t18, t28, t31 * t16 - t28 * t32 + 0; t28 * t19, -t28 * t18, -t31, t28 * t16 + t31 * t32 + 0; t18, t19, 0, t38; 0, 0, 0, 1; t31 * t13, -t10, t28, t33; t28 * t13, -t9, -t31, t37; t12, t13, 0, t36; 0, 0, 0, 1; t4, t3, t10, t35 * t31 + t33; t2, t1, t9, t35 * t28 + t37; t8, -t44, -t13, t12 * pkin(4) - t13 * pkin(9) + t36; 0, 0, 0, 1; t4, t3, t10, t34 * t31 + (pkin(5) * t26 - t22) * t28 + t39; t2, t1, t9, -pkin(5) * t41 + t34 * t28 + t37; t8, -t44, -t13, t12 * t15 + t13 * t25 + t36; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
