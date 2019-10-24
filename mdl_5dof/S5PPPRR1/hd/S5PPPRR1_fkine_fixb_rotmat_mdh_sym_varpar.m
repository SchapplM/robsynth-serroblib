% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:17:12
% EndTime: 2019-10-24 10:17:12
% DurationCPUTime: 0.15s
% Computational Cost: add. (128->52), mult. (153->63), div. (0->0), fcn. (216->10), ass. (0->33)
t20 = pkin(9) + qJ(4);
t14 = sin(t20);
t22 = sin(pkin(8));
t41 = t22 * t14;
t15 = cos(t20);
t40 = t22 * t15;
t21 = sin(pkin(9));
t23 = sin(pkin(7));
t39 = t23 * t21;
t10 = t23 * t22;
t25 = cos(pkin(8));
t38 = t23 * t25;
t26 = cos(pkin(7));
t11 = t26 * t22;
t37 = t26 * t25;
t36 = t23 * pkin(1) + 0;
t19 = qJ(1) + 0;
t35 = t26 * pkin(1) + t23 * qJ(2) + 0;
t24 = cos(pkin(9));
t13 = t24 * pkin(3) + pkin(2);
t27 = -pkin(5) - qJ(3);
t34 = t22 * t13 + t25 * t27 + t19;
t33 = pkin(2) * t25 + qJ(3) * t22;
t32 = -t26 * qJ(2) + t36;
t31 = pkin(3) * t39 - t27 * t11 + t13 * t37 + t35;
t30 = -t27 * t10 + t13 * t38 + (-pkin(3) * t21 - qJ(2)) * t26 + t36;
t29 = cos(qJ(5));
t28 = sin(qJ(5));
t4 = t23 * t14 + t15 * t37;
t3 = t14 * t37 - t23 * t15;
t2 = -t26 * t14 + t15 * t38;
t1 = t14 * t38 + t26 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t37, -t11, t23, t35; t38, -t10, -t26, t32; t22, t25, 0, t19; 0, 0, 0, 1; t24 * t37 + t39, -t21 * t37 + t23 * t24, t11, t33 * t26 + t35; -t26 * t21 + t24 * t38, -t21 * t38 - t26 * t24, t10, t33 * t23 + t32; t22 * t24, -t22 * t21, -t25, t22 * pkin(2) - t25 * qJ(3) + t19; 0, 0, 0, 1; t4, -t3, t11, t31; t2, -t1, t10, t30; t40, -t41, -t25, t34; 0, 0, 0, 1; t28 * t11 + t4 * t29, t29 * t11 - t4 * t28, t3, t4 * pkin(4) + t3 * pkin(6) + t31; t28 * t10 + t2 * t29, t29 * t10 - t2 * t28, t1, t2 * pkin(4) + t1 * pkin(6) + t30; -t25 * t28 + t29 * t40, -t25 * t29 - t28 * t40, t41, (pkin(4) * t15 + pkin(6) * t14) * t22 + t34; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
