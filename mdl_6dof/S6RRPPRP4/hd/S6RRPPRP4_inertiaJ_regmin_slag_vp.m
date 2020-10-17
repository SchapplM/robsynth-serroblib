% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:21:49
% EndTime: 2019-05-06 09:21:52
% DurationCPUTime: 0.78s
% Computational Cost: add. (695->118), mult. (1347->212), div. (0->0), fcn. (1398->6), ass. (0->65)
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t81 = t47 ^ 2 + t48 ^ 2;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t27 = t51 * t47 - t49 * t48;
t50 = sin(qJ(2));
t21 = t27 * t50;
t80 = 0.2e1 * t21;
t26 = t49 * t47 + t51 * t48;
t79 = 0.2e1 * t26;
t78 = -0.2e1 * t27;
t56 = t47 * qJ(4) + pkin(2);
t29 = -t48 * pkin(3) - t56;
t77 = -0.2e1 * t29;
t76 = 0.2e1 * t50;
t52 = cos(qJ(2));
t75 = 0.2e1 * t52;
t74 = pkin(3) + pkin(4);
t73 = pkin(2) * t47;
t72 = pkin(7) * t48;
t71 = t50 * pkin(7);
t70 = t52 * pkin(5);
t69 = t52 * pkin(7);
t30 = -t52 * pkin(2) - t50 * qJ(3) - pkin(1);
t35 = t47 * t69;
t43 = t52 * pkin(3);
t10 = t52 * pkin(4) + t35 + t43 + (-pkin(8) * t50 - t30) * t48;
t19 = t47 * t30 + t48 * t69;
t16 = -t52 * qJ(4) + t19;
t37 = t47 * t50;
t12 = pkin(8) * t37 + t16;
t4 = t49 * t10 + t51 * t12;
t63 = -pkin(8) + qJ(3);
t31 = t63 * t48;
t57 = t63 * t47;
t13 = t49 * t31 - t51 * t57;
t68 = t13 * t52;
t14 = t51 * t31 + t49 * t57;
t67 = t14 * t52;
t39 = t48 * t50;
t65 = t49 * t52;
t62 = t81 * qJ(3) ^ 2;
t61 = qJ(3) * t52;
t60 = t47 * qJ(3);
t59 = t52 * qJ(6);
t58 = -t51 * t10 + t49 * t12;
t18 = t48 * t30 - t35;
t17 = -t18 + t43;
t55 = t16 * t48 + t17 * t47;
t54 = -t18 * t47 + t19 * t48;
t23 = t48 * t74 + t56;
t34 = qJ(4) * t39;
t15 = t34 + (-t47 * t74 - pkin(7)) * t50;
t46 = t50 ^ 2;
t40 = t51 * t52;
t33 = t52 * t60;
t28 = 0.2e1 * t81 * qJ(3);
t22 = t26 * t50;
t20 = -t34 + (pkin(3) * t47 + pkin(7)) * t50;
t7 = t26 * pkin(5) - t27 * qJ(6) + t23;
t5 = t21 * pkin(5) + t22 * qJ(6) - t15;
t2 = t58 - t70;
t1 = t59 + t4;
t3 = [1, 0, 0, t46, t50 * t75, 0, 0, 0, pkin(1) * t75, -0.2e1 * pkin(1) * t50, 0.2e1 * t46 * pkin(7) * t47 - 0.2e1 * t18 * t52, 0.2e1 * t19 * t52 + 0.2e1 * t46 * t72 (-t18 * t48 - t19 * t47) * t76, t46 * pkin(7) ^ 2 + t18 ^ 2 + t19 ^ 2, 0.2e1 * t17 * t52 + 0.2e1 * t20 * t37 (-t16 * t47 + t17 * t48) * t76, -0.2e1 * t16 * t52 - 0.2e1 * t20 * t39, t16 ^ 2 + t17 ^ 2 + t20 ^ 2, t22 ^ 2, t22 * t80, t22 * t75, t52 * t80, t52 ^ 2, -0.2e1 * t15 * t21 - 0.2e1 * t52 * t58, 0.2e1 * t15 * t22 - 0.2e1 * t4 * t52, -0.2e1 * t2 * t52 + 0.2e1 * t5 * t21, 0.2e1 * t1 * t21 + 0.2e1 * t2 * t22, 0.2e1 * t1 * t52 + 0.2e1 * t5 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t50, t52, 0, -t71, -t69, t33 + (-t72 - t73) * t50, pkin(7) * t37 + (-pkin(2) * t50 + t61) * t48, t54, -pkin(2) * t71 + qJ(3) * t54, -t20 * t48 + t29 * t37 + t33, t55, -t20 * t47 + (-t29 * t50 - t61) * t48, qJ(3) * t55 + t20 * t29, t22 * t27, t27 * t21 - t22 * t26, t27 * t52, -t26 * t52, 0, t15 * t26 - t23 * t21 - t68, t15 * t27 + t23 * t22 - t67, -t7 * t21 - t5 * t26 - t68, -t1 * t26 + t13 * t22 + t14 * t21 + t2 * t27, -t7 * t22 + t5 * t27 + t67, t1 * t14 + t2 * t13 - t5 * t7; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t48, -0.2e1 * t73, t28, pkin(2) ^ 2 + t62, t48 * t77, t28, t47 * t77, t29 ^ 2 + t62, t27 ^ 2, t26 * t78, 0, 0, 0, t23 * t79, 0.2e1 * t23 * t27, t7 * t79, 0.2e1 * t13 * t27 - 0.2e1 * t14 * t26, t7 * t78, t13 ^ 2 + t14 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t71, t37, 0, -t39, t20, 0, 0, 0, 0, 0, t21, -t22, t21, 0, t22, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, -pkin(2), -t48, 0, -t47, t29, 0, 0, 0, 0, 0, -t26, -t27, -t26, 0, t27, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t39, 0, t17, 0, 0, 0, 0, 0, t40, -t65, t40, t49 * t21 - t51 * t22, t65, t1 * t49 - t2 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, -t49 * t26 - t51 * t27, 0, -t13 * t51 + t14 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 ^ 2 + t51 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t21, t52, -t58, -t4, -t58 + 0.2e1 * t70, -t22 * pkin(5) + t21 * qJ(6), 0.2e1 * t59 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, -t13, -t14, -t13, -pkin(5) * t27 - t26 * qJ(6), t14, -t13 * pkin(5) + t14 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, t51, 0, t49, t51 * pkin(5) + t49 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
