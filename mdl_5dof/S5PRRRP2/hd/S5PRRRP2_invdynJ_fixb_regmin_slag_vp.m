% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:42:00
% EndTime: 2019-12-05 16:42:04
% DurationCPUTime: 0.77s
% Computational Cost: add. (926->173), mult. (1293->194), div. (0->0), fcn. (662->8), ass. (0->112)
t64 = pkin(8) + qJ(2);
t61 = qJ(3) + t64;
t50 = cos(t61);
t138 = g(2) * t50;
t49 = sin(t61);
t45 = g(1) * t49;
t146 = t138 - t45;
t68 = sin(qJ(4));
t70 = cos(qJ(4));
t86 = t70 * pkin(4) + t68 * qJ(5);
t31 = -pkin(3) - t86;
t69 = sin(qJ(3));
t107 = qJDD(2) * t69;
t71 = cos(qJ(3));
t115 = qJD(3) * t71;
t63 = qJDD(2) + qJDD(3);
t145 = qJD(1) * qJD(4) + t63 * pkin(7) + (qJD(2) * t115 + t107) * pkin(2);
t65 = qJD(2) + qJD(3);
t120 = pkin(2) * qJD(2);
t99 = t69 * t120;
t28 = t65 * pkin(7) + t99;
t127 = t68 * t28;
t16 = t70 * qJD(1) - t127;
t144 = qJD(5) - t16;
t142 = g(1) * t50 + g(2) * t49;
t102 = t68 * qJDD(1) + t145 * t70;
t104 = qJDD(4) * qJ(5);
t4 = t104 + (qJD(5) - t127) * qJD(4) + t102;
t111 = qJDD(4) * pkin(4);
t140 = qJDD(5) - t111;
t112 = qJD(4) * t70;
t93 = -t70 * qJDD(1) + t28 * t112 + t145 * t68;
t5 = t93 + t140;
t143 = t4 * t70 + t5 * t68;
t10 = -qJD(4) * pkin(4) + t144;
t66 = t68 ^ 2;
t67 = t70 ^ 2;
t121 = t66 + t67;
t141 = t121 * t65;
t108 = qJD(4) * qJ(5);
t17 = t68 * qJD(1) + t70 * t28;
t11 = t17 + t108;
t72 = qJD(4) ^ 2;
t139 = pkin(7) * t72;
t137 = t63 * pkin(3);
t136 = t65 * pkin(3);
t134 = t71 * pkin(2);
t133 = t31 * t63;
t132 = t31 * t65;
t131 = t49 * t68;
t130 = t50 * t68;
t53 = t69 * pkin(2) + pkin(7);
t129 = t53 * t72;
t128 = t65 * t68;
t126 = t70 * t63;
t125 = g(1) * t131 - g(2) * t130;
t123 = -qJD(3) * t99 + qJDD(2) * t134;
t122 = t66 - t67;
t118 = pkin(7) * qJDD(4);
t117 = qJD(2) * t71;
t116 = qJD(3) * t69;
t114 = qJD(4) * t65;
t113 = qJD(4) * t68;
t110 = t11 * qJD(4);
t109 = t68 * qJD(5);
t106 = qJDD(4) * t53;
t62 = t65 ^ 2;
t103 = t68 * t62 * t70;
t38 = t70 * t45;
t91 = t65 * t99;
t98 = pkin(2) * t117;
t101 = t98 * t113 + t70 * t91 + t38;
t100 = pkin(2) * t115;
t97 = t65 * t116;
t18 = -t123 - t137;
t96 = -t18 - t138;
t95 = t121 * t63;
t94 = t16 + t127;
t29 = -t98 - t136;
t92 = t29 * t112 + t18 * t68 - t125;
t89 = -t137 + t139;
t57 = sin(t64);
t58 = cos(t64);
t87 = g(1) * t57 - g(2) * t58;
t85 = pkin(4) * t68 - qJ(5) * t70;
t84 = t10 * t68 + t11 * t70;
t83 = t123 - t146;
t1 = (t85 * qJD(4) - t109) * t65 + t133 - t123;
t82 = -t1 - t133 - t139;
t26 = t31 - t134;
t81 = t26 * t65 - t100;
t21 = pkin(4) * t113 - t70 * t108 - t109;
t80 = t10 * t112 - t68 * t110 - t142 + t143;
t79 = g(1) * t130 + g(2) * t131 - g(3) * t70 - t93;
t12 = pkin(2) * t116 + t21;
t78 = -t12 * t65 - t26 * t63 - t1 - t129;
t54 = -pkin(3) - t134;
t77 = pkin(2) * t97 + t54 * t63 + t129;
t76 = t17 * qJD(4) + t79;
t75 = -t106 + (t54 * t65 - t100) * qJD(4);
t74 = (t10 * t70 - t11 * t68) * qJD(4) + t143;
t73 = -t142 * pkin(7) + t146 * t31;
t47 = t68 * t63;
t33 = qJDD(4) * t70 - t72 * t68;
t32 = qJDD(4) * t68 + t72 * t70;
t23 = t29 * t113;
t22 = t85 * t65;
t20 = 0.2e1 * t112 * t128 + t66 * t63;
t9 = -t98 + t132;
t7 = -0.2e1 * t122 * t114 + 0.2e1 * t68 * t126;
t6 = t9 * t113;
t2 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t33, 0, t32, t84 * qJD(4) + t4 * t68 - t5 * t70 - g(3); 0, qJDD(2), t87, g(1) * t58 + g(2) * t57, t63, (t63 * t71 - t97) * pkin(2) + t83, ((-qJDD(2) - t63) * t69 + (-qJD(2) - t65) * t115) * pkin(2) + t142, t20, t7, t32, t33, 0, t23 + t38 + t75 * t68 + (-t77 + t96) * t70, t77 * t68 + t75 * t70 + t92, t38 + t6 + (t81 * qJD(4) - t106) * t68 + (t78 - t138) * t70, t100 * t141 + t53 * t95 + t80, (t106 + (-t81 - t9) * qJD(4)) * t70 + t78 * t68 + t125, t1 * t26 + t9 * t12 + (t84 * t115 + t87) * pkin(2) + t74 * t53 + t73; 0, 0, 0, 0, t63, t83 + t91, (-t107 + (-qJD(3) + t65) * t117) * pkin(2) + t142, t20, t7, t32, t33, 0, t23 + (-pkin(3) * t114 - t118) * t68 + (-t89 + t96) * t70 + t101, (-t118 + (t98 - t136) * qJD(4)) * t70 + (t89 - t91) * t68 + t92, t6 + (t31 * t114 - t118) * t68 + (-t21 * t65 - t138 + t82) * t70 + t101, pkin(7) * t95 - t98 * t141 + t80, (t118 + (-t9 - t98 - t132) * qJD(4)) * t70 + ((-t21 + t99) * t65 + t82) * t68 + t125, t1 * t31 + t9 * t21 + (-t69 * t9 - t84 * t71) * t120 + t74 * pkin(7) + t73; 0, 0, 0, 0, 0, 0, 0, -t103, t122 * t62, t47, t126, qJDD(4), -t29 * t128 + t76, g(3) * t68 + t94 * qJD(4) + (-t29 * t65 + t142) * t70 - t102, 0.2e1 * t111 - qJDD(5) + (t22 * t70 - t68 * t9) * t65 + t76, -t85 * t63, 0.2e1 * t104 + (t22 * t65 - g(3)) * t68 + (t65 * t9 - t142) * t70 + (0.2e1 * qJD(5) - t94) * qJD(4) + t102, -t5 * pkin(4) - g(3) * t86 + t4 * qJ(5) - t10 * t17 + t144 * t11 + t142 * t85 - t9 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t103, t47, -t66 * t62 - t72, t9 * t128 - t110 + t140 - t79;];
tau_reg = t2;
