% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(pkin(10));
t47 = cos(pkin(10));
t50 = sin(qJ(6));
t52 = cos(qJ(6));
t28 = t52 * t45 + t50 * t47;
t46 = sin(pkin(9));
t48 = cos(pkin(9));
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t29 = t51 * t46 - t53 * t48;
t16 = t28 * t29;
t75 = 0.2e1 * t16;
t30 = -t50 * t45 + t52 * t47;
t31 = t53 * t46 + t51 * t48;
t74 = -0.2e1 * t31;
t37 = t45 * pkin(5) + qJ(4);
t73 = 0.2e1 * t37;
t39 = -t48 * pkin(2) - pkin(1);
t72 = 0.2e1 * t39;
t71 = 0.2e1 * qJ(4);
t63 = pkin(3) + qJ(5);
t70 = -pkin(8) - t63;
t56 = -t31 * qJ(4) + t39;
t10 = t63 * t29 + t56;
t62 = pkin(7) + qJ(2);
t34 = t62 * t46;
t35 = t62 * t48;
t20 = t53 * t34 + t51 * t35;
t13 = t31 * pkin(4) + t20;
t7 = t47 * t10 + t45 * t13;
t23 = t28 * t31;
t24 = t30 * t31;
t69 = t45 * t29;
t68 = t45 * t31;
t67 = t47 * t29;
t66 = t47 * t31;
t61 = t45 ^ 2 + t47 ^ 2;
t60 = t46 ^ 2 + t48 ^ 2;
t59 = qJ(4) * t29;
t12 = t47 * t13;
t6 = -t45 * t10 + t12;
t3 = t7 * t45 + t6 * t47;
t58 = -t6 * t45 + t7 * t47;
t21 = -t51 * t34 + t53 * t35;
t57 = t31 * t63 + t59;
t54 = qJ(4) ^ 2;
t33 = t70 * t47;
t32 = t70 * t45;
t27 = t61 * t63;
t26 = t31 ^ 2;
t19 = t52 * t32 + t50 * t33;
t18 = -t50 * t32 + t52 * t33;
t17 = t29 * pkin(3) + t56;
t15 = t30 * t29;
t14 = -t29 * pkin(4) + t21;
t8 = (-pkin(5) * t47 - pkin(4)) * t29 + t21;
t5 = pkin(8) * t67 + t7;
t4 = t31 * pkin(5) + t12 + (-pkin(8) * t29 - t10) * t45;
t2 = t50 * t4 + t52 * t5;
t1 = t52 * t4 - t50 * t5;
t9 = [1, 0, 0, 0.2e1 * pkin(1) * t48, -0.2e1 * pkin(1) * t46, 0.2e1 * t60 * qJ(2), t60 * qJ(2) ^ 2 + pkin(1) ^ 2, t26, t29 * t74, 0, 0, 0, t29 * t72, t31 * t72, 0.2e1 * t20 * t31 - 0.2e1 * t21 * t29, -0.2e1 * t17 * t29, t17 * t74, t17 ^ 2 + t20 ^ 2 + t21 ^ 2, -0.2e1 * t14 * t67 + 0.2e1 * t6 * t31, 0.2e1 * t14 * t69 - 0.2e1 * t7 * t31, 0.2e1 * t58 * t29, t14 ^ 2 + t6 ^ 2 + t7 ^ 2, t16 ^ 2, t15 * t75, t31 * t75, -t15 * t74, t26, 0.2e1 * t1 * t31 - 0.2e1 * t8 * t15, 0.2e1 * t8 * t16 - 0.2e1 * t2 * t31; 0, 0, 0, -t48, t46, 0, -pkin(1), 0, 0, 0, 0, 0, t29, t31, 0, -t29, -t31, t17, -t68, -t66, t61 * t29, t58, 0, 0, 0, 0, 0, -t23, -t24; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, -t20, -t21, -pkin(3) * t31 - t59, t20, t21, -t20 * pkin(3) + t21 * qJ(4), t14 * t45 - t57 * t47, t14 * t47 + t57 * t45, -t3, t14 * qJ(4) - t3 * t63, t16 * t30, t30 * t15 - t16 * t28, t24, -t23, 0, -t37 * t15 + t18 * t31 + t8 * t28, t37 * t16 - t19 * t31 + t8 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t71, pkin(3) ^ 2 + t54, t45 * t71, t47 * t71, 0.2e1 * t27, t61 * t63 ^ 2 + t54, t30 ^ 2, -0.2e1 * t30 * t28, 0, 0, 0, t28 * t73, t30 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, t20, t66, -t68, 0, t3, 0, 0, 0, 0, 0, t24, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -t61, -t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t69, 0, t14, 0, 0, 0, 0, 0, -t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t47, 0, qJ(4), 0, 0, 0, 0, 0, t28, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
