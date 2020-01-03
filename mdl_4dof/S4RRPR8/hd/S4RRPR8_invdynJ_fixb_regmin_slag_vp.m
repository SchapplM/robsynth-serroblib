% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:24
% EndTime: 2019-12-31 17:08:27
% DurationCPUTime: 1.10s
% Computational Cost: add. (566->178), mult. (1273->238), div. (0->0), fcn. (802->6), ass. (0->120)
t65 = sin(qJ(4));
t69 = cos(qJ(2));
t126 = t65 * t69;
t66 = sin(qJ(2));
t68 = cos(qJ(4));
t145 = t66 * t68 - t126;
t144 = qJD(2) * t126 + qJD(4) * t145;
t109 = qJD(1) * qJD(2);
t104 = t66 * t109;
t110 = t69 * qJDD(1);
t141 = t104 - t110;
t103 = t69 * t109;
t55 = t66 * qJDD(1);
t142 = t103 + t55;
t26 = t65 * t66 + t68 * t69;
t81 = t26 * qJD(4);
t1 = -qJD(1) * t81 + t141 * t65 + t142 * t68;
t20 = t26 * qJD(1);
t59 = qJD(2) - qJD(4);
t143 = -t20 * t59 + t1;
t2 = qJD(1) * t144 + t26 * qJDD(1) - t68 * t104;
t115 = qJD(1) * t69;
t116 = qJD(1) * t66;
t22 = -t115 * t65 + t116 * t68;
t140 = t22 * t59 + t2;
t138 = qJD(4) + t59;
t70 = cos(qJ(1));
t130 = g(2) * t70;
t67 = sin(qJ(1));
t133 = g(1) * t67;
t137 = -t130 + t133;
t131 = g(2) * t67;
t97 = g(1) * t70 + t131;
t118 = pkin(5) * qJDD(2);
t56 = t66 * qJ(3);
t102 = pkin(1) + t56;
t57 = t69 * pkin(2);
t86 = t102 + t57;
t24 = t86 * qJD(1);
t121 = t57 + t56;
t34 = -pkin(1) - t121;
t136 = (qJD(1) * t34 - t24) * qJD(2) - t118;
t135 = pkin(2) + pkin(3);
t134 = pkin(5) - pkin(6);
t128 = t22 * t20;
t124 = t66 * t70;
t73 = qJD(1) ^ 2;
t123 = t66 * t73;
t63 = t66 ^ 2;
t64 = t69 ^ 2;
t120 = t63 - t64;
t119 = qJ(3) * t69;
t112 = t66 * qJD(3);
t96 = pkin(2) * t66 - t119;
t19 = qJD(2) * t96 - t112;
t117 = qJD(1) * t19;
t114 = qJD(2) * t66;
t113 = qJDD(2) * pkin(2);
t53 = pkin(5) * t116;
t111 = -pkin(6) * t116 + qJD(3) + t53;
t62 = qJD(2) * qJ(3);
t108 = t69 * t123;
t37 = t134 * t69;
t106 = t20 ^ 2 - t22 ^ 2;
t47 = pkin(5) * t103;
t51 = pkin(5) * t55;
t105 = qJDD(3) + t47 + t51;
t100 = -qJD(2) * pkin(2) + qJD(3);
t99 = t59 ^ 2;
t54 = pkin(5) * t115;
t31 = -pkin(6) * t115 + t54;
t72 = qJD(2) ^ 2;
t98 = pkin(5) * t72 + t130;
t12 = -qJD(2) * t135 + t111;
t23 = t31 + t62;
t95 = t68 * t12 - t65 * t23;
t94 = -t65 * t12 - t68 * t23;
t33 = t100 + t53;
t35 = t54 + t62;
t93 = t33 * t69 - t35 * t66;
t36 = t134 * t66;
t92 = t36 * t68 - t37 * t65;
t91 = t36 * t65 + t37 * t68;
t89 = qJ(3) * t68 - t135 * t65;
t88 = -qJ(3) * t65 - t135 * t68;
t87 = g(1) * t124 - g(3) * t69 + t131 * t66 - t51;
t85 = -t135 * t66 + t119;
t84 = -qJDD(3) + t87;
t83 = -0.2e1 * pkin(1) * t109 - t118;
t80 = t135 * t69 + t102;
t78 = 0.2e1 * qJDD(1) * pkin(1) - t98;
t52 = pkin(5) * t110;
t60 = qJDD(2) * qJ(3);
t61 = qJD(2) * qJD(3);
t11 = -pkin(5) * t104 + t52 + t60 + t61;
t4 = -qJDD(1) * t86 + t117;
t77 = -qJDD(1) * t34 - t117 - t4 - t98;
t9 = qJD(2) * t85 + t112;
t10 = t80 * qJD(1);
t15 = t145 * t67;
t17 = -t124 * t68 + t126 * t70;
t5 = -pkin(6) * t142 - t135 * qJDD(2) + t105;
t6 = pkin(6) * t141 + t11;
t76 = g(1) * t17 - g(2) * t15 + g(3) * t26 - t10 * t22 + t68 * t5 - t65 * t6;
t16 = t26 * t67;
t18 = t26 * t70;
t75 = g(1) * t18 + g(2) * t16 + g(3) * t145 + t10 * t20 - t65 * t5 - t68 * t6;
t13 = t105 - t113;
t74 = qJD(2) * t93 + t11 * t69 + t13 * t66 - t97;
t58 = qJDD(2) - qJDD(4);
t49 = t69 * t133;
t32 = qJD(2) * t37;
t30 = t134 * t114;
t28 = t96 * qJD(1);
t25 = pkin(3) * t69 - t34;
t14 = t85 * qJD(1);
t8 = qJD(2) * t26 - t81;
t7 = -t114 * t68 + t144;
t3 = qJD(1) * t9 + qJDD(1) * t80;
t21 = [qJDD(1), t137, t97, qJDD(1) * t63 + 0.2e1 * t103 * t66, -0.2e1 * t109 * t120 + 0.2e1 * t110 * t66, qJDD(2) * t66 + t69 * t72, qJDD(2) * t69 - t66 * t72, 0, t66 * t83 + t69 * t78 + t49, t83 * t69 + (-t78 - t133) * t66, t136 * t66 + t69 * t77 + t49, (t63 + t64) * qJDD(1) * pkin(5) + t74, -t136 * t69 + (t77 + t133) * t66, t74 * pkin(5) + t137 * t86 - t24 * t19 + t4 * t34, t1 * t145 + t22 * t8, -t1 * t26 - t145 * t2 - t20 * t8 - t22 * t7, -t145 * t58 - t59 * t8, t26 * t58 + t59 * t7, 0, t9 * t20 + t25 * t2 + t3 * t26 + t10 * t7 - (-qJD(4) * t91 + t65 * t30 + t68 * t32) * t59 - t92 * t58 + g(1) * t16 - g(2) * t18, t9 * t22 + t25 * t1 + t3 * t145 + t10 * t8 + (qJD(4) * t92 - t68 * t30 + t65 * t32) * t59 + t91 * t58 + g(1) * t15 + g(2) * t17; 0, 0, 0, -t108, t120 * t73, t55, t110, qJDD(2), pkin(1) * t123 + t87, g(3) * t66 - t52 + (pkin(1) * t73 + t97) * t69, 0.2e1 * t113 + (t24 * t66 + t28 * t69) * qJD(1) + t84, -t96 * qJDD(1) + ((t35 - t62) * t66 + (t100 - t33) * t69) * qJD(1), t52 + 0.2e1 * t60 + 0.2e1 * t61 + (qJD(1) * t28 - g(3)) * t66 + (-qJD(1) * t24 - t97) * t69, -t93 * qJD(1) * pkin(5) - t13 * pkin(2) - g(3) * t121 + t11 * qJ(3) + t35 * qJD(3) + t24 * t28 + t96 * t97, -t128, t106, -t143, t140, t58, -t88 * t58 - t14 * t20 + (t111 * t65 + t68 * t31) * t59 + (t59 * t89 - t94) * qJD(4) - t76, t89 * t58 - t14 * t22 + (t111 * t68 - t65 * t31) * t59 + (t59 * t88 + t95) * qJD(4) - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t108, t55, -t63 * t73 - t72, -qJD(2) * t35 - t116 * t24 - t113 + t47 - t84, 0, 0, 0, 0, 0, -t116 * t20 - t68 * t58 - t65 * t99, -t116 * t22 + t65 * t58 - t68 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t106, t143, -t140, -t58, t138 * t94 + t76, -t138 * t95 + t75;];
tau_reg = t21;
