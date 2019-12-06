% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:17
% EndTime: 2019-12-05 16:52:21
% DurationCPUTime: 0.83s
% Computational Cost: add. (1220->164), mult. (2999->234), div. (0->0), fcn. (2022->6), ass. (0->109)
t76 = sin(qJ(4));
t77 = sin(qJ(3));
t79 = cos(qJ(4));
t80 = cos(qJ(3));
t53 = t76 * t80 + t79 * t77;
t73 = qJD(3) + qJD(4);
t142 = t73 * t53;
t21 = t142 * qJD(2);
t111 = qJD(2) * qJD(3);
t141 = -0.2e1 * t111;
t128 = t79 * t80;
t130 = t76 * t77;
t52 = -t128 + t130;
t78 = sin(qJ(2));
t45 = t52 * t78;
t82 = qJD(3) ^ 2;
t83 = qJD(2) ^ 2;
t140 = (t82 + t83) * t78;
t137 = -pkin(7) - pkin(6);
t105 = qJD(3) * t137;
t54 = t77 * t105;
t55 = t80 * t105;
t81 = cos(qJ(2));
t87 = t52 * t81;
t56 = t137 * t77;
t57 = t137 * t80;
t96 = t79 * t56 + t76 * t57;
t135 = -qJD(1) * t87 - t96 * qJD(4) - t79 * t54 - t76 * t55;
t139 = t135 * t73;
t114 = t78 * qJD(1);
t118 = qJD(3) * t77;
t89 = pkin(3) * t118 - t114;
t48 = t53 * qJD(2);
t138 = t48 ^ 2;
t35 = t76 * t56 - t79 * t57;
t88 = t53 * t81;
t136 = qJD(1) * t88 - t35 * qJD(4) - t76 * t54 + t79 * t55;
t60 = qJD(2) * pkin(6) + t114;
t103 = pkin(7) * qJD(2) + t60;
t43 = t103 * t80;
t131 = t76 * t43;
t42 = t103 * t77;
t38 = qJD(3) * pkin(3) - t42;
t15 = t79 * t38 - t131;
t134 = t15 * t73;
t129 = t79 * t43;
t16 = t76 * t38 + t129;
t133 = t16 * t73;
t107 = qJD(2) * t128;
t120 = qJD(2) * t77;
t108 = t76 * t120;
t46 = -t107 + t108;
t132 = t48 * t46;
t127 = t82 * t77;
t126 = t82 * t80;
t115 = qJD(4) * t79;
t18 = -t79 * t42 - t131;
t125 = pkin(3) * t115 + qJD(5) - t18;
t104 = t80 * t111;
t124 = -qJD(4) * t107 - t79 * t104;
t109 = pkin(3) * t120;
t51 = qJD(2) * t114 + qJD(3) * t109;
t123 = t77 ^ 2 - t80 ^ 2;
t121 = qJD(2) * pkin(2);
t119 = qJD(2) * t78;
t117 = qJD(3) * t80;
t116 = qJD(4) * t76;
t113 = t81 * qJD(1);
t112 = qJD(5) - t15;
t71 = -t80 * pkin(3) - pkin(2);
t106 = t136 * t73;
t97 = t73 * t130;
t28 = -t80 * t115 - t79 * t117 + t97;
t102 = -pkin(4) * t142 - t28 * qJ(5) + t53 * qJD(5) - t89;
t30 = -t60 * t118 + (-pkin(7) * t118 + t80 * t113) * qJD(2);
t31 = -t60 * t117 + (-pkin(7) * t117 - t77 * t113) * qJD(2);
t101 = -t38 * t115 + t43 * t116 - t79 * t30 - t76 * t31;
t2 = t43 * t115 + t38 * t116 + t76 * t30 - t79 * t31;
t99 = t81 * t141;
t17 = -t76 * t42 + t129;
t98 = pkin(3) * t116 - t17;
t22 = t48 * pkin(4) + t46 * qJ(5);
t95 = qJD(2) * t121;
t6 = qJD(2) * t88 - t73 * t45;
t94 = t46 * t119 - t81 * t21 - t6 * t73;
t50 = t71 * qJD(2) - t113;
t13 = t46 * pkin(4) - t48 * qJ(5) + t50;
t93 = -t13 * t48 - t2;
t92 = -t13 * t46 - t101;
t91 = -t50 * t48 - t2;
t90 = t50 * t46 + t101;
t86 = -0.2e1 * qJD(3) * t121;
t20 = qJD(2) * t97 + t124;
t5 = -qJD(2) * t87 - t142 * t78;
t85 = -t48 * t119 - t81 * t20 + t5 * t73;
t72 = t73 * qJD(5);
t70 = -t79 * pkin(3) - pkin(4);
t67 = t76 * pkin(3) + qJ(5);
t44 = t53 * t78;
t23 = t52 * pkin(4) - t53 * qJ(5) + t71;
t19 = t22 + t109;
t14 = -t46 ^ 2 + t138;
t12 = t73 * qJ(5) + t16;
t11 = -t73 * pkin(4) + t112;
t10 = t48 * t73 - t21;
t9 = -t124 + (-t108 + t46) * t73;
t3 = t21 * pkin(4) + t20 * qJ(5) - t48 * qJD(5) + t51;
t1 = t72 - t101;
t4 = [0, 0, -t83 * t78, -t83 * t81, 0, 0, 0, 0, 0, -t80 * t140 + t77 * t99, t77 * t140 + t80 * t99, 0, 0, 0, 0, 0, t94, -t85, t94, -t44 * t20 + t45 * t21 - t5 * t46 + t6 * t48, t85, -t1 * t45 + t11 * t6 + t13 * t119 + t12 * t5 + t2 * t44 - t3 * t81; 0, 0, 0, 0, 0.2e1 * t77 * t104, t123 * t141, t126, -t127, 0, -pkin(6) * t126 + t77 * t86, pkin(6) * t127 + t80 * t86, -t20 * t53 - t48 * t28, -t142 * t48 + t20 * t52 - t53 * t21 + t28 * t46, -t28 * t73, -t142 * t73, 0, t142 * t50 + t71 * t21 + t89 * t46 + t51 * t52 + t106, -t71 * t20 - t50 * t28 + t89 * t48 + t51 * t53 + t139, -t102 * t46 + t13 * t142 + t23 * t21 + t3 * t52 + t106, -t1 * t52 - t11 * t28 - t12 * t142 + t135 * t46 - t136 * t48 + t2 * t53 + t20 * t96 - t35 * t21, t102 * t48 + t13 * t28 + t23 * t20 - t3 * t53 - t139, t1 * t35 - t102 * t13 - t136 * t11 - t135 * t12 - t2 * t96 + t3 * t23; 0, 0, 0, 0, -t77 * t83 * t80, t123 * t83, 0, 0, 0, t77 * t95, t80 * t95, t132, t14, t9, t10, 0, t17 * t73 + (-t73 * t116 - t46 * t120) * pkin(3) + t91, t18 * t73 + (-t73 * t115 - t48 * t120) * pkin(3) + t90, -t19 * t46 - t98 * t73 + t93, -t70 * t20 - t67 * t21 + (t12 + t98) * t48 + (t11 - t125) * t46, t125 * t73 + t19 * t48 + t72 + t92, t1 * t67 + t98 * t11 + t125 * t12 - t13 * t19 + t2 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t14, t9, t10, 0, t91 + t133, t90 + t134, -t22 * t46 + t133 + t93, pkin(4) * t20 - t21 * qJ(5) + (t12 - t16) * t48 + (t11 - t112) * t46, t22 * t48 - t134 + 0.2e1 * t72 + t92, -t2 * pkin(4) + t1 * qJ(5) - t11 * t16 + t112 * t12 - t13 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t9, -t73 ^ 2 - t138, -t12 * t73 - t93;];
tauc_reg = t4;
