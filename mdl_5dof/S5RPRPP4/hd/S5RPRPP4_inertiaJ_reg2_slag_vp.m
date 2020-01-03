% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = sin(pkin(7));
t29 = cos(pkin(7));
t30 = sin(qJ(3));
t31 = cos(qJ(3));
t12 = t28 * t31 + t29 * t30;
t14 = -t28 * t30 + t29 * t31;
t50 = (t12 * t28 + t14 * t29) * pkin(3);
t10 = t12 ^ 2;
t11 = t14 ^ 2;
t49 = t10 + t11;
t32 = -pkin(1) - pkin(6);
t41 = -qJ(4) + t32;
t17 = t41 * t30;
t38 = t41 * t31;
t6 = t28 * t17 - t29 * t38;
t8 = t29 * t17 + t28 * t38;
t39 = t8 * t12 - t6 * t14;
t23 = t30 * pkin(3) + qJ(2);
t48 = 0.2e1 * t23;
t47 = 0.2e1 * qJ(2);
t46 = t28 * pkin(3);
t45 = t29 * pkin(3);
t42 = t14 * t12;
t25 = t30 ^ 2;
t26 = t31 ^ 2;
t18 = t25 + t26;
t40 = t6 ^ 2 + t8 ^ 2;
t19 = qJ(5) + t46;
t21 = pkin(4) + t45;
t37 = t12 * t19 + t21 * t14;
t35 = 0.2e1 * t39;
t33 = qJ(2) ^ 2;
t16 = t18 * t32;
t4 = t12 * pkin(4) - t14 * qJ(5) + t23;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t47, (pkin(1) ^ 2) + t33, t26, -0.2e1 * t31 * t30, 0, t25, 0, 0, t30 * t47, t31 * t47, -0.2e1 * t16, t18 * t32 ^ 2 + t33, t11, -0.2e1 * t42, 0, t10, 0, 0, t12 * t48, t14 * t48, -t35, t23 ^ 2 + t40, t11, 0, 0.2e1 * t42, 0, 0, t10, 0.2e1 * t4 * t12, -t35, -0.2e1 * t4 * t14, t4 ^ 2 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t18, t16, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t39, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t30, 0, t31 * t32, -t30 * t32, 0, 0, 0, 0, t14, 0, -t12, 0, -t6, -t8, -t50, (t28 * t8 - t29 * t6) * pkin(3), 0, t14, 0, 0, t12, 0, -t6, -t37, t8, t8 * t19 - t6 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, 0, t50, 0, 0, 0, 0, 0, 0, t14, 0, t12, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t46, 0, (t28 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t21, 0, 0.2e1 * t19, t19 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14, 0, t23, 0, 0, 0, 0, 0, 0, t12, 0, -t14, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
