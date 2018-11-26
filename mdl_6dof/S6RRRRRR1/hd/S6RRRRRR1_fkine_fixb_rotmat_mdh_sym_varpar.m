% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:37:12
% EndTime: 2018-11-23 18:37:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (198->53), mult. (89->54), div. (0->0), fcn. (138->12), ass. (0->35)
t30 = -pkin(8) - pkin(7);
t26 = sin(qJ(1));
t23 = qJ(2) + qJ(3);
t17 = qJ(4) + t23;
t14 = qJ(5) + t17;
t7 = sin(t14);
t42 = t26 * t7;
t29 = cos(qJ(1));
t41 = t29 * t7;
t28 = cos(qJ(2));
t13 = t28 * pkin(2) + pkin(1);
t24 = sin(qJ(6));
t40 = t26 * t24;
t27 = cos(qJ(6));
t39 = t26 * t27;
t38 = t29 * t24;
t37 = t29 * t27;
t22 = -pkin(9) + t30;
t21 = pkin(6) + 0;
t16 = cos(t23);
t4 = pkin(3) * t16 + t13;
t18 = -pkin(10) + t22;
t12 = cos(t17);
t3 = pkin(4) * t12 + t4;
t36 = t29 * t18 + t26 * t3 + 0;
t25 = sin(qJ(2));
t35 = t25 * pkin(2) + t21;
t15 = sin(t23);
t34 = pkin(3) * t15 + t35;
t8 = cos(t14);
t33 = pkin(5) * t8 + pkin(11) * t7;
t11 = sin(t17);
t32 = pkin(4) * t11 + t34;
t31 = -t26 * t18 + t29 * t3 + 0;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t29 * t28, -t29 * t25, t26, t29 * pkin(1) + t26 * pkin(7) + 0; t26 * t28, -t26 * t25, -t29, t26 * pkin(1) - t29 * pkin(7) + 0; t25, t28, 0, t21; 0, 0, 0, 1; t29 * t16, -t29 * t15, t26, t29 * t13 - t26 * t30 + 0; t26 * t16, -t26 * t15, -t29, t26 * t13 + t29 * t30 + 0; t15, t16, 0, t35; 0, 0, 0, 1; t29 * t12, -t29 * t11, t26, -t26 * t22 + t29 * t4 + 0; t26 * t12, -t26 * t11, -t29, t29 * t22 + t26 * t4 + 0; t11, t12, 0, t34; 0, 0, 0, 1; t29 * t8, -t41, t26, t31; t26 * t8, -t42, -t29, t36; t7, t8, 0, t32; 0, 0, 0, 1; t8 * t37 + t40, -t8 * t38 + t39, t41, t33 * t29 + t31; t8 * t39 - t38, -t8 * t40 - t37, t42, t33 * t26 + t36; t7 * t27, -t7 * t24, -t8, t7 * pkin(5) - t8 * pkin(11) + t32; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
