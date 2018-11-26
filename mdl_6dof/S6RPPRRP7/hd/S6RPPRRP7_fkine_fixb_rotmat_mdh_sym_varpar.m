% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2018-11-23 15:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRRP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:47
% EndTime: 2018-11-23 15:46:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (131->51), mult. (110->46), div. (0->0), fcn. (160->8), ass. (0->39)
t15 = pkin(9) + qJ(4);
t10 = cos(t15);
t21 = sin(qJ(5));
t43 = t10 * t21;
t22 = sin(qJ(1));
t42 = t22 * t10;
t17 = sin(pkin(9));
t41 = t22 * t17;
t40 = t22 * t21;
t23 = cos(qJ(5));
t39 = t22 * t23;
t24 = cos(qJ(1));
t38 = t24 * t21;
t37 = t24 * t23;
t16 = pkin(6) + 0;
t36 = t22 * pkin(1) + 0;
t35 = pkin(2) + t16;
t20 = -pkin(7) - qJ(3);
t34 = pkin(5) * t21 - t20;
t33 = -pkin(3) * t17 - qJ(2);
t32 = t24 * pkin(1) + t22 * qJ(2) + 0;
t18 = cos(pkin(9));
t31 = t18 * pkin(3) + t35;
t30 = pkin(3) * t41 + t32;
t9 = sin(t15);
t29 = pkin(4) * t9 - pkin(8) * t10;
t19 = -qJ(6) - pkin(8);
t8 = t23 * pkin(5) + pkin(4);
t28 = t10 * t19 + t8 * t9;
t27 = -t22 * t20 + t36;
t26 = -t24 * qJ(2) + t36;
t25 = -t24 * t20 + t30;
t6 = t24 * t10;
t5 = t10 * t23;
t4 = -t9 * t37 + t40;
t3 = t9 * t38 + t39;
t2 = t9 * t39 + t38;
t1 = -t9 * t40 + t37;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t22, 0, 0; t22, t24, 0, 0; 0, 0, 1, t16; 0, 0, 0, 1; 0, -t24, t22, t32; 0, -t22, -t24, t26; 1, 0, 0, t16; 0, 0, 0, 1; t41, t22 * t18, t24, t24 * qJ(3) + t32; -t24 * t17, -t24 * t18, t22, t22 * qJ(3) + t26; t18, -t17, 0, t35; 0, 0, 0, 1; t22 * t9, t42, t24, t25; -t24 * t9, -t6, t22, t33 * t24 + t27; t10, -t9, 0, t31; 0, 0, 0, 1; t2, t1, -t42, t29 * t22 + t25; t4, t3, t6 (-t29 + t33) * t24 + t27; t5, -t43, t9, t10 * pkin(4) + t9 * pkin(8) + t31; 0, 0, 0, 1; t2, t1, -t42, t28 * t22 + t34 * t24 + t30; t4, t3, t6, t34 * t22 + (-t28 + t33) * t24 + t36; t5, -t43, t9, t10 * t8 - t9 * t19 + t31; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
