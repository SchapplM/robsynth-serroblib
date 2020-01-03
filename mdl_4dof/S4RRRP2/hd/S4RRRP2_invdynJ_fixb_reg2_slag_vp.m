% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:12
% EndTime: 2019-12-31 17:13:13
% DurationCPUTime: 0.91s
% Computational Cost: add. (930->200), mult. (1449->254), div. (0->0), fcn. (720->8), ass. (0->134)
t75 = sin(qJ(3));
t71 = t75 ^ 2;
t78 = cos(qJ(3));
t72 = t78 ^ 2;
t137 = t71 + t72;
t76 = sin(qJ(2));
t126 = qJDD(1) * t76;
t79 = cos(qJ(2));
t131 = qJD(2) * t79;
t69 = qJDD(1) + qJDD(2);
t19 = t69 * pkin(6) + (qJD(1) * t131 + t126) * pkin(1);
t111 = t137 * t19;
t157 = t137 * t79;
t136 = pkin(1) * qJD(1);
t119 = t76 * t136;
t70 = qJD(1) + qJD(2);
t31 = t70 * pkin(6) + t119;
t133 = qJD(1) * t79;
t118 = pkin(1) * t133;
t32 = -t70 * pkin(2) - t118;
t158 = t157 * t31 + t32 * t76;
t73 = qJ(1) + qJ(2);
t63 = sin(t73);
t55 = g(2) * t63;
t64 = cos(t73);
t57 = g(1) * t64;
t140 = t55 + t57;
t105 = qJ(4) * t70 + t31;
t12 = t105 * t75;
t13 = t105 * t78;
t129 = qJD(3) * t75;
t116 = t70 * t129;
t59 = t78 * pkin(3) + pkin(2);
t147 = t59 * t69;
t156 = -pkin(3) * t116 + t147;
t155 = pkin(3) * t71;
t56 = g(1) * t63;
t77 = sin(qJ(1));
t154 = g(1) * t77;
t153 = g(2) * t64;
t152 = g(3) * t78;
t151 = t69 * pkin(2);
t150 = t79 * pkin(1);
t135 = qJD(3) * pkin(3);
t9 = -t12 + t135;
t149 = t12 + t9;
t146 = t63 * t78;
t145 = t64 * t75;
t68 = t70 ^ 2;
t144 = t68 * t78;
t143 = t70 * t75;
t142 = t70 * t78;
t49 = t75 * t69;
t50 = t78 * t69;
t74 = -qJ(4) - pkin(6);
t141 = t64 * pkin(2) + t63 * pkin(6);
t139 = -qJD(2) * t119 + qJDD(1) * t150;
t138 = t71 - t72;
t58 = t76 * pkin(1) + pkin(6);
t134 = -qJ(4) - t58;
t132 = qJD(2) * t76;
t130 = qJD(3) * t70;
t128 = qJD(3) * t78;
t127 = qJDD(3) * pkin(3);
t106 = -qJ(4) * t69 - t19;
t87 = qJD(4) * t70 - t106;
t90 = qJD(3) * t105;
t3 = -t75 * t90 + t78 * t87;
t125 = t3 * t78 - t140;
t17 = -t59 * t70 + qJD(4) - t118;
t45 = g(2) * t145;
t5 = qJDD(4) - t139 - t156;
t124 = t17 * t128 + t5 * t75 + t45;
t18 = -t139 - t151;
t123 = t32 * t128 + t18 * t75 + t45;
t101 = t70 * t119;
t46 = g(1) * t146;
t99 = qJD(3) * t118;
t122 = t78 * t101 + t75 * t99 + t46;
t121 = g(2) * t146 + g(3) * t75 + t78 * t57;
t120 = pkin(1) * t131;
t117 = t70 * t132;
t10 = t17 * t129;
t115 = t70 * t128;
t114 = -t5 - t153;
t113 = -t18 - t153;
t110 = t137 * t69;
t108 = -t32 * t70 - t19;
t107 = t64 * t59 - t63 * t74;
t104 = qJD(3) * t74;
t103 = -t140 + t111;
t102 = qJD(3) * t134;
t100 = t75 * t115;
t98 = g(1) * t145 + t75 * t55 - t152;
t97 = g(1) * (-t63 * pkin(2) + t64 * pkin(6));
t81 = qJD(3) ^ 2;
t96 = -pkin(6) * t81 + t151;
t80 = cos(qJ(1));
t95 = -g(2) * t80 + t154;
t94 = -t13 * t78 + t75 * t9;
t29 = pkin(1) * t132 + pkin(3) * t129;
t36 = -t59 - t150;
t93 = t29 * t70 + t36 * t69;
t92 = -t63 * t59 - t64 * t74;
t91 = -t139 - t56 + t153;
t89 = t157 * t136;
t88 = -pkin(2) * t130 - pkin(6) * qJDD(3);
t86 = -t101 - t56;
t85 = (-qJD(4) - t17) * t70 + t106;
t60 = -pkin(2) - t150;
t84 = pkin(1) * t117 + t58 * t81 + t60 * t69;
t83 = -qJDD(3) * t58 + (t60 * t70 - t120) * qJD(3);
t66 = t80 * pkin(1);
t65 = t78 * qJ(4);
t62 = t78 * qJD(4);
t41 = t75 * t144;
t40 = t78 * pkin(6) + t65;
t39 = t74 * t75;
t38 = t78 * t99;
t35 = qJDD(3) * t78 - t81 * t75;
t34 = qJDD(3) * t75 + t81 * t78;
t28 = t78 * t58 + t65;
t27 = t138 * t68;
t26 = t134 * t75;
t24 = t32 * t129;
t23 = -t75 * qJD(4) + t78 * t104;
t22 = t75 * t104 + t62;
t21 = t72 * t69 - 0.2e1 * t100;
t20 = t71 * t69 + 0.2e1 * t100;
t8 = -0.2e1 * t138 * t130 + 0.2e1 * t75 * t50;
t7 = (-qJD(4) - t120) * t75 + t78 * t102;
t6 = t75 * t102 + t78 * t120 + t62;
t2 = -t75 * t87 - t78 * t90 + t127;
t1 = [0, 0, 0, 0, 0, qJDD(1), t95, g(1) * t80 + g(2) * t77, 0, 0, 0, 0, 0, 0, 0, t69, (t69 * t79 - t117) * pkin(1) - t91, ((-qJDD(1) - t69) * t76 + (-qJD(1) - t70) * t131) * pkin(1) + t140, 0, (t95 + (t76 ^ 2 + t79 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t20, t8, t34, t21, t35, 0, t24 + t46 + t83 * t75 + (t113 - t84) * t78, t83 * t78 + (t84 - t56) * t75 + t123, pkin(1) * qJD(2) * t157 * t70 + t58 * t110 + t103, t18 * t60 - t97 - g(2) * (t66 + t141) + t58 * t111 + (t158 * qJD(2) + t154) * pkin(1), t20, t8, t34, t21, t35, 0, t26 * qJDD(3) + t10 + t46 + (t36 * t143 + t7) * qJD(3) + (t114 - t93) * t78, -t28 * qJDD(3) + (t36 * t142 - t6) * qJD(3) + (t93 - t56) * t75 + t124, (t28 * t69 + t6 * t70 + (-t26 * t70 - t9) * qJD(3)) * t78 + (-t26 * t69 - t7 * t70 - t2 + (-t28 * t70 - t13) * qJD(3)) * t75 + t125, t3 * t28 + t13 * t6 + t2 * t26 + t9 * t7 + t5 * t36 + t17 * t29 - g(1) * (-t77 * pkin(1) + t92) - g(2) * (t107 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t91 + t101, (-t126 + (-qJD(2) + t70) * t133) * pkin(1) + t140, 0, 0, t20, t8, t34, t21, t35, 0, t24 + t88 * t75 + (t113 + t96) * t78 + t122, t38 + t88 * t78 + (t86 - t96) * t75 + t123, pkin(6) * t110 - t70 * t89 + t103, -t18 * pkin(2) + pkin(6) * t111 - g(2) * t141 - t158 * t136 - t97, t20, t8, t34, t21, t35, 0, t39 * qJDD(3) + t10 + (-t59 * t143 + t23) * qJD(3) + (t114 + t156) * t78 + t122, -t40 * qJDD(3) + t38 + (t86 - t147) * t75 + (-t22 + (-t59 * t78 + t155) * t70) * qJD(3) + t124, (-qJD(3) * t9 + t40 * t69) * t78 + (-t13 * qJD(3) - t39 * t69 - t2) * t75 + (t22 * t78 - t23 * t75 + (-t39 * t78 - t40 * t75) * qJD(3) - t89) * t70 + t125, t3 * t40 + t13 * t22 + t2 * t39 + t9 * t23 - t5 * t59 + pkin(3) * t10 - g(1) * t92 - g(2) * t107 + (-t17 * t76 + t79 * t94) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t27, t49, t41, t50, qJDD(3), t108 * t75 + t98, t108 * t78 + t121, 0, 0, -t41, t27, t49, t41, t50, qJDD(3), 0.2e1 * t127 + (pkin(3) * t144 + t85) * t75 + t98, -t68 * t155 + t85 * t78 + t121, -pkin(3) * t49 + (-t135 + t149) * t142, t149 * t13 + (-t152 + t2 + (-t17 * t70 + t140) * t75) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 + 0.2e1 * t116, t49 + 0.2e1 * t115, -t137 * t68, t70 * t94 - t114 - t56;];
tau_reg = t1;
