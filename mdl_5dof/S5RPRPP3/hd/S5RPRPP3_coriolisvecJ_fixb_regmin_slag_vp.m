% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:53
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.59s
% Computational Cost: add. (986->171), mult. (2630->193), div. (0->0), fcn. (1817->4), ass. (0->97)
t120 = pkin(6) + qJ(2);
t80 = sin(pkin(7));
t64 = t120 * t80;
t61 = qJD(1) * t64;
t81 = cos(pkin(7));
t65 = t120 * t81;
t62 = qJD(1) * t65;
t83 = sin(qJ(3));
t84 = cos(qJ(3));
t119 = -t84 * t61 - t83 * t62;
t110 = -qJD(4) + t119;
t60 = t84 * t80 + t83 * t81;
t56 = t60 * qJD(3);
t37 = qJD(1) * t56;
t122 = t84 * t81;
t105 = qJD(1) * t122;
t123 = t83 * t80;
t106 = qJD(1) * t123;
t51 = -t105 + t106;
t129 = t37 * qJ(5) + t51 * qJD(5);
t127 = t51 ^ 2;
t53 = t60 * qJD(1);
t46 = t53 ^ 2;
t128 = -t127 - t46;
t113 = qJD(3) * t84;
t19 = (qJD(2) * t80 + qJD(3) * t65) * t83 - qJD(2) * t122 + t64 * t113;
t126 = t51 * pkin(4);
t75 = -t81 * pkin(2) - pkin(1);
t63 = t75 * qJD(1) + qJD(2);
t88 = -t53 * qJ(4) + t63;
t18 = t51 * pkin(3) + t88;
t125 = t18 * t53;
t124 = t53 * t51;
t121 = pkin(3) + qJ(5);
t31 = -t83 * t61 + t84 * t62;
t117 = t80 ^ 2 + t81 ^ 2;
t116 = qJ(4) * t37;
t115 = t51 * qJ(4);
t114 = qJD(3) * t83;
t112 = t19 * qJD(3);
t99 = -t83 * t64 + t84 * t65;
t20 = t60 * qJD(2) + t99 * qJD(3);
t111 = t20 * qJD(3);
t16 = -t53 * pkin(4) + t119;
t109 = qJD(4) - t16;
t17 = t31 - t126;
t108 = -qJD(5) - t17;
t107 = qJD(1) * qJD(2);
t104 = t117 * qJD(1) ^ 2;
t103 = t83 * t107;
t102 = t84 * t107;
t32 = t84 * t64 + t83 * t65;
t101 = t81 * t102 - t80 * t103 - t61 * t113 - t62 * t114;
t13 = t80 * t102 + t81 * t103 + t62 * t113 - t61 * t114;
t27 = -qJD(3) * qJ(4) - t31;
t70 = qJD(3) * t105;
t36 = qJD(3) * t106 - t70;
t100 = t37 * pkin(3) + t36 * qJ(4);
t55 = -t81 * t113 + t80 * t114;
t98 = t55 * qJ(4) - t60 * qJD(4);
t97 = 0.2e1 * t117 * t107;
t95 = -t60 * qJ(4) + t75;
t94 = t36 * pkin(4) - t13;
t93 = -t37 * pkin(4) + t101;
t92 = qJD(3) * t119 - t101;
t91 = t31 * qJD(3) - t13;
t8 = -t53 * qJD(4) + t100;
t5 = t121 * t51 + t88;
t90 = t5 * t53 - t94;
t89 = -t5 * t51 + t93;
t26 = 0.2e1 * t53 * qJD(3);
t85 = qJD(3) ^ 2;
t79 = qJD(3) * qJD(4);
t76 = 0.2e1 * t79;
t59 = -t122 + t123;
t39 = qJD(3) * t51;
t35 = -t85 - t46;
t29 = t59 * pkin(3) + t95;
t28 = t53 * pkin(3) + t115;
t25 = t70 + (t51 - t106) * qJD(3);
t24 = -t70 + (t51 + t106) * qJD(3);
t23 = -qJD(3) * pkin(3) - t110;
t22 = -t59 * pkin(4) + t99;
t21 = t60 * pkin(4) + t32;
t15 = t121 * t59 + t95;
t14 = t56 * pkin(3) + t98;
t12 = t121 * t53 + t115;
t11 = qJD(5) - t27 - t126;
t10 = -t79 - t101;
t9 = -t121 * qJD(3) + t109;
t7 = -t55 * pkin(4) + t20;
t6 = -t56 * pkin(4) - t19;
t4 = -qJD(3) * qJD(5) - t94;
t3 = t79 + t93;
t2 = t59 * qJD(5) + t121 * t56 + t98;
t1 = t8 + t129;
t30 = [0, 0, 0, 0, 0, t97, qJ(2) * t97, -t36 * t60 - t53 * t55, t36 * t59 - t60 * t37 + t55 * t51 - t53 * t56, -t55 * qJD(3), -t56 * qJD(3), 0, t75 * t37 + t63 * t56 - t111, -t75 * t36 - t63 * t55 + t112, t10 * t59 + t13 * t60 + t19 * t51 + t20 * t53 - t23 * t55 + t27 * t56 - t32 * t36 - t37 * t99, -t14 * t51 - t18 * t56 - t29 * t37 - t8 * t59 + t111, -t14 * t53 + t18 * t55 + t29 * t36 - t8 * t60 - t112, -t10 * t99 + t13 * t32 + t18 * t14 + t27 * t19 + t23 * t20 + t8 * t29, -t11 * t56 - t21 * t36 - t22 * t37 - t3 * t59 + t4 * t60 - t6 * t51 + t7 * t53 - t9 * t55, t6 * qJD(3) - t1 * t60 + t15 * t36 - t2 * t53 + t5 * t55, -t7 * qJD(3) + t1 * t59 + t15 * t37 + t2 * t51 + t5 * t56, t1 * t15 + t11 * t6 + t5 * t2 + t4 * t21 + t3 * t22 + t9 * t7; 0, 0, 0, 0, 0, -t104, -qJ(2) * t104, 0, 0, 0, 0, 0, t26, -t24, t128, -t26, t36 + t39, -t27 * t51 + (-qJD(4) - t23) * t53 + t100, t128, t24, t26, t11 * t51 + (-qJD(4) - t9) * t53 + t100 + t129; 0, 0, 0, 0, 0, 0, 0, t124, -t127 + t46, t25, 0, 0, -t63 * t53 + t91, t63 * t51 + t92, pkin(3) * t36 - t116 + (-t27 - t31) * t53 + (t23 + t110) * t51, t28 * t51 + t125 - t91, -t18 * t51 + t28 * t53 + t76 - t92, -t13 * pkin(3) - t10 * qJ(4) + t110 * t27 - t18 * t28 - t23 * t31, -t116 + t121 * t36 + (t11 + t108) * t53 + (t9 - t109) * t51, -t16 * qJD(3) + t12 * t53 + t76 + t89, -t12 * t51 + (0.2e1 * qJD(5) + t17) * qJD(3) - t90, t3 * qJ(4) + t108 * t9 + t109 * t11 - t5 * t12 - t121 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 + t39, -t124, t35, t27 * qJD(3) + t125 + t13, t25, t35, t124, (-qJD(5) - t11) * qJD(3) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t127 - t85, t9 * qJD(3) + t79 + t89;];
tauc_reg = t30;
