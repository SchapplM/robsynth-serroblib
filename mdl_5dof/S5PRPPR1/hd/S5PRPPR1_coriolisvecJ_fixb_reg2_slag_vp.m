% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:10
% EndTime: 2019-12-05 15:22:13
% DurationCPUTime: 0.62s
% Computational Cost: add. (947->125), mult. (2514->202), div. (0->0), fcn. (1835->6), ass. (0->84)
t66 = sin(qJ(5));
t62 = sin(pkin(9));
t63 = sin(pkin(8));
t90 = qJD(2) * t63;
t79 = t62 * t90;
t64 = cos(pkin(9));
t67 = cos(qJ(5));
t97 = t67 * t64;
t81 = t63 * t97;
t26 = qJD(2) * t81 - t66 * t79;
t65 = cos(pkin(8));
t57 = t65 * qJD(2) - qJD(5);
t87 = qJD(4) * t63;
t88 = qJD(3) * t65;
t38 = -t62 * t88 - t64 * t87;
t33 = qJD(2) * t38;
t39 = -t62 * t87 + t64 * t88;
t34 = qJD(2) * t39;
t52 = -t65 * pkin(3) - t63 * qJ(4) - pkin(2);
t37 = qJD(2) * t52 + qJD(3);
t85 = qJ(3) * qJD(2);
t50 = t63 * qJD(1) + t65 * t85;
t14 = t62 * t37 + t64 * t50;
t10 = -pkin(6) * t79 + t14;
t13 = t64 * t37 - t62 * t50;
t105 = pkin(6) * t63;
t83 = t64 * t105;
t9 = (-t65 * pkin(4) - t83) * qJD(2) + t13;
t75 = t66 * t10 - t67 * t9;
t1 = -qJD(5) * t75 + t66 * t33 + t67 * t34;
t106 = t26 ^ 2;
t47 = t67 * t62 + t66 * t64;
t35 = t47 * t63;
t29 = qJD(5) * t35;
t18 = qJD(2) * t29;
t104 = t18 * t65;
t70 = qJD(2) * t47;
t23 = t63 * t70;
t103 = t23 * t57;
t102 = t26 * t23;
t101 = t26 * t57;
t100 = t29 * t57;
t19 = t26 * qJD(5);
t99 = t65 * t19;
t98 = t66 * t62;
t46 = t97 - t98;
t36 = t46 * t63;
t96 = -t36 * t19 + t29 * t23;
t95 = t47 * qJD(5) - t65 * t70;
t94 = t57 * t46;
t91 = qJ(3) * t65;
t93 = t62 * t52 + t64 * t91;
t60 = t63 ^ 2;
t61 = t65 ^ 2;
t92 = t60 + t61;
t89 = qJD(3) * t63;
t84 = qJD(2) * qJD(3);
t68 = qJD(2) ^ 2;
t82 = t63 * t65 * t68;
t80 = 0.2e1 * qJD(3) * t60;
t78 = t92 * t68;
t58 = t63 * t84;
t49 = t65 * qJD(1) - t63 * t85;
t6 = t67 * t10 + t66 * t9;
t44 = t64 * t52;
t15 = -t83 + t44 + (-qJ(3) * t62 - pkin(4)) * t65;
t16 = -t62 * t105 + t93;
t7 = t67 * t15 - t66 * t16;
t8 = t66 * t15 + t67 * t16;
t30 = (-t63 * t98 + t81) * qJD(5);
t74 = -t35 * t18 + t26 * t30;
t73 = -t33 * t64 - t34 * t62;
t71 = t49 * t63 - t50 * t65;
t45 = qJD(4) - t49;
t69 = qJD(5) * t23;
t2 = -qJD(5) * t6 + t67 * t33 - t66 * t34;
t56 = t60 * qJ(3) * t84;
t48 = (pkin(4) * t62 + qJ(3)) * t63;
t22 = pkin(4) * t79 + t45;
t21 = t23 ^ 2;
t17 = t30 * t57;
t4 = -qJD(5) * t8 + t67 * t38 - t66 * t39;
t3 = qJD(5) * t7 + t66 * t38 + t67 * t39;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t33 * t62 + t34 * t64 - t65 * t84) * t63, 0, 0, 0, 0, 0, 0, t17 - t99, -t100 + t104, t74 + t96, t1 * t36 - t2 * t35 - t6 * t29 + t30 * t75 - t65 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t92 * t84, t56 + (t61 * t85 - t71) * qJD(3), 0, 0, 0, 0, 0, 0, -t33 * t65 + (-t38 * t65 + t62 * t80) * qJD(2), t34 * t65 + (t39 * t65 + t64 * t80) * qJD(2), ((-t38 * t64 - t39 * t62) * qJD(2) + t73) * t63, t34 * t93 + t14 * t39 + t33 * (-t62 * t91 + t44) + t13 * t38 + t56 + t45 * t89, -t18 * t36 - t26 * t29, -t74 + t96, t100 + t104, t19 * t35 + t23 * t30, t17 + t99, 0, t48 * t19 - t2 * t65 + t22 * t30 - t4 * t57 + (qJD(2) * t35 + t23) * t89, t1 * t65 - t48 * t18 - t22 * t29 + t3 * t57 + (qJD(2) * t36 + t26) * t89, -t1 * t35 + t7 * t18 - t8 * t19 - t2 * t36 - t3 * t23 - t4 * t26 - t29 * t75 - t6 * t30, t1 * t8 + t2 * t7 + t6 * t3 - t75 * t4 + (qJD(2) * t48 + t22) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t71 * qJD(2), 0, 0, 0, 0, 0, 0, -t62 * t78, -t64 * t78, 0, (-t45 * t63 + (t13 * t62 - t14 * t64) * t65) * qJD(2) - t73, 0, 0, 0, 0, 0, 0, -t23 * t90 + t95 * t57, -t26 * t90 - t94 * t57, t46 * t18 - t47 * t19 + t94 * t23 + t95 * t26, t1 * t47 + t2 * t46 - t22 * t90 - t94 * t6 + t75 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t82, t62 * t82, (-t62 ^ 2 - t64 ^ 2) * t68 * t60, t58 + (t13 * t64 + t14 * t62) * t90, 0, 0, 0, 0, 0, 0, t19 - t101, -t69 + t103, -t21 - t106, t6 * t23 - t26 * t75 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t21 + t106, -t69 - t103, -t102, -t19 - t101, 0, -t22 * t26 - t6 * t57 + t2, t22 * t23 + t57 * t75 - t1, 0, 0;];
tauc_reg = t5;
