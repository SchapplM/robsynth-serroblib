% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:50:35
% EndTime: 2018-11-23 17:50:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (198->53), mult. (89->54), div. (0->0), fcn. (138->12), ass. (0->35)
t30 = -pkin(8) - pkin(7);
t26 = sin(qJ(1));
t23 = qJ(2) + qJ(3);
t15 = pkin(11) + t23;
t13 = qJ(5) + t15;
t7 = sin(t13);
t42 = t26 * t7;
t29 = cos(qJ(1));
t41 = t29 * t7;
t28 = cos(qJ(2));
t14 = t28 * pkin(2) + pkin(1);
t24 = sin(qJ(6));
t40 = t26 * t24;
t27 = cos(qJ(6));
t39 = t26 * t27;
t38 = t29 * t24;
t37 = t29 * t27;
t22 = pkin(6) + 0;
t21 = -qJ(4) + t30;
t17 = cos(t23);
t4 = pkin(3) * t17 + t14;
t18 = -pkin(9) + t21;
t10 = cos(t15);
t3 = pkin(4) * t10 + t4;
t36 = t29 * t18 + t26 * t3 + 0;
t25 = sin(qJ(2));
t35 = t25 * pkin(2) + t22;
t16 = sin(t23);
t34 = pkin(3) * t16 + t35;
t8 = cos(t13);
t33 = pkin(5) * t8 + pkin(10) * t7;
t32 = -t26 * t18 + t29 * t3 + 0;
t9 = sin(t15);
t31 = pkin(4) * t9 + t34;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t29 * t28, -t29 * t25, t26, t29 * pkin(1) + t26 * pkin(7) + 0; t26 * t28, -t26 * t25, -t29, t26 * pkin(1) - t29 * pkin(7) + 0; t25, t28, 0, t22; 0, 0, 0, 1; t29 * t17, -t29 * t16, t26, t29 * t14 - t26 * t30 + 0; t26 * t17, -t26 * t16, -t29, t26 * t14 + t29 * t30 + 0; t16, t17, 0, t35; 0, 0, 0, 1; t29 * t10, -t29 * t9, t26, -t26 * t21 + t29 * t4 + 0; t26 * t10, -t26 * t9, -t29, t29 * t21 + t26 * t4 + 0; t9, t10, 0, t34; 0, 0, 0, 1; t29 * t8, -t41, t26, t32; t26 * t8, -t42, -t29, t36; t7, t8, 0, t31; 0, 0, 0, 1; t8 * t37 + t40, -t8 * t38 + t39, t41, t33 * t29 + t32; t8 * t39 - t38, -t8 * t40 - t37, t42, t33 * t26 + t36; t7 * t27, -t7 * t24, -t8, t7 * pkin(5) - t8 * pkin(10) + t31; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
