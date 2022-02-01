% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:55
% DurationCPUTime: 0.76s
% Computational Cost: add. (1247->166), mult. (2126->223), div. (0->0), fcn. (1241->16), ass. (0->121)
t92 = cos(qJ(2));
t131 = qJD(2) * t92;
t119 = qJD(1) * t131;
t88 = sin(qJ(2));
t124 = qJDD(1) * t88;
t101 = (t119 + t124) * pkin(1);
t133 = pkin(1) * qJD(1);
t122 = t88 * t133;
t149 = t92 * pkin(1);
t72 = qJDD(1) * t149;
t79 = qJDD(1) + qJDD(2);
t36 = t79 * pkin(2) - qJD(2) * t122 + t72;
t80 = qJD(1) + qJD(2);
t48 = t80 * pkin(2) + t92 * t133;
t128 = qJD(3) * t88;
t118 = qJD(1) * t128;
t114 = pkin(1) * t118;
t87 = sin(qJ(3));
t51 = t87 * t114;
t91 = cos(qJ(3));
t15 = (qJD(3) * t48 + t101) * t91 + t87 * t36 - t51;
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t74 = qJDD(3) + t79;
t129 = qJD(3) * t87;
t32 = t91 * t36;
t111 = -t48 * t129 + t32;
t127 = qJD(3) * t91;
t121 = t88 * t127;
t95 = (-t87 * t124 + (-t87 * t131 - t121) * qJD(1)) * pkin(1) + t111;
t9 = t74 * pkin(3) + t95;
t4 = -t84 * t15 + t85 * t9;
t83 = qJ(1) + qJ(2);
t78 = qJ(3) + t83;
t67 = pkin(9) + t78;
t55 = cos(t67);
t158 = -t74 * pkin(4) + g(2) * t55 - t4;
t68 = sin(t78);
t69 = cos(t78);
t156 = g(1) * t69 + g(2) * t68;
t155 = g(1) * t68 - g(2) * t69;
t132 = pkin(2) * qJD(3);
t145 = t85 * t87;
t143 = t88 * t91;
t108 = -t87 * t92 - t143;
t42 = t108 * t133;
t144 = t87 * t88;
t107 = t91 * t92 - t144;
t43 = t107 * t133;
t139 = -t85 * t42 + t84 * t43 - (t84 * t91 + t145) * t132;
t147 = t84 * t87;
t70 = t91 * pkin(2) + pkin(3);
t106 = -pkin(2) * t147 + t85 * t70;
t37 = -pkin(4) - t106;
t136 = pkin(2) * t145 + t84 * t70;
t38 = pkin(8) + t136;
t75 = qJD(3) + t80;
t94 = qJD(5) ^ 2;
t154 = t139 * t75 - t37 * t74 - t38 * t94;
t153 = pkin(2) * t74;
t54 = sin(t67);
t152 = g(1) * t54;
t5 = t85 * t15 + t84 * t9;
t31 = t91 * t122 + t87 * t48;
t148 = t84 * t31;
t146 = t85 * t31;
t90 = cos(qJ(5));
t142 = t90 * t74;
t30 = -t87 * t122 + t91 * t48;
t28 = t75 * pkin(3) + t30;
t16 = t85 * t28 - t148;
t13 = -t75 * pkin(4) - t16;
t86 = sin(qJ(5));
t141 = t13 * qJD(5) * t86 + t90 * t152;
t71 = pkin(2) + t149;
t53 = t91 * t71;
t41 = -pkin(1) * t144 + pkin(3) + t53;
t44 = pkin(1) * t143 + t87 * t71;
t140 = t84 * t41 + t85 * t44;
t138 = -t84 * t42 - t85 * t43 + (t85 * t91 - t147) * t132;
t77 = cos(t83);
t137 = pkin(2) * t77 + pkin(3) * t69;
t76 = sin(t83);
t135 = g(1) * t77 + g(2) * t76;
t81 = t86 ^ 2;
t134 = -t90 ^ 2 + t81;
t126 = t90 * qJD(5);
t125 = qJDD(4) - g(3);
t123 = t13 * t126 + t158 * t86;
t116 = qJD(1) * (-qJD(2) + t80);
t115 = qJD(2) * (-qJD(1) - t80);
t113 = g(1) * t76 - g(2) * t77 + t72;
t112 = -pkin(2) * t76 - pkin(3) * t68;
t109 = t85 * t41 - t84 * t44;
t20 = -pkin(4) - t109;
t21 = pkin(8) + t140;
t24 = t71 * t127 + (t107 * qJD(2) - t87 * t128) * pkin(1);
t25 = -t71 * t129 + (t108 * qJD(2) - t121) * pkin(1);
t6 = t84 * t24 - t85 * t25;
t104 = t20 * t74 + t21 * t94 + t6 * t75;
t18 = t84 * t30 + t146;
t61 = t84 * pkin(3) + pkin(8);
t62 = -t85 * pkin(3) - pkin(4);
t103 = -t18 * t75 + t61 * t94 + t62 * t74;
t102 = -t74 * pkin(8) + g(1) * t55 + g(2) * t54 - t13 * t75 - t5;
t7 = t85 * t24 + t84 * t25;
t100 = -qJDD(5) * t21 + (t20 * t75 - t7) * qJD(5);
t19 = t85 * t30 - t148;
t99 = -qJDD(5) * t61 + (t62 * t75 + t19) * qJD(5);
t98 = -qJDD(5) * t38 + (t37 * t75 - t138) * qJD(5);
t97 = (-pkin(2) * t75 - t48) * qJD(3) - t101;
t96 = -t15 + t156;
t93 = cos(qJ(1));
t89 = sin(qJ(1));
t73 = t75 ^ 2;
t50 = qJDD(5) * t90 - t94 * t86;
t49 = qJDD(5) * t86 + t94 * t90;
t33 = 0.2e1 * t86 * t75 * t126 + t81 * t74;
t26 = -0.2e1 * t134 * t75 * qJD(5) + 0.2e1 * t86 * t142;
t17 = t84 * t28 + t146;
t1 = [qJDD(1), g(1) * t89 - g(2) * t93, g(1) * t93 + g(2) * t89, t79, (t88 * t115 + t79 * t92) * pkin(1) + t113, ((-qJDD(1) - t79) * t88 + t92 * t115) * pkin(1) + t135, t74, t25 * t75 + t53 * t74 + (-t91 * t118 + (-t119 + (-qJDD(1) - t74) * t88) * t87) * pkin(1) + t111 + t155, -t24 * t75 - t44 * t74 + t96, t5 * t140 + t17 * t7 + t4 * t109 - t16 * t6 - g(1) * (-t89 * pkin(1) + t112) - g(2) * (t93 * pkin(1) + t137), t33, t26, t49, t50, 0, t100 * t86 + (-t104 - t158) * t90 + t141, t100 * t90 + (t104 - t152) * t86 + t123; 0, 0, 0, t79, t88 * pkin(1) * t116 + t113, (t92 * t116 - t124) * pkin(1) + t135, t74, -t42 * t75 + t32 + (-t114 + t153) * t91 + t97 * t87 + t155, t43 * t75 + t51 + (-t36 - t153) * t87 + t97 * t91 + t156, -g(1) * t112 - g(2) * t137 + t4 * t106 + t5 * t136 + t138 * t17 + t139 * t16, t33, t26, t49, t50, 0, t98 * t86 + (-t158 + t154) * t90 + t141, t98 * t90 + (-t152 - t154) * t86 + t123; 0, 0, 0, 0, 0, 0, t74, t31 * t75 + t155 + t95, t30 * t75 + t96, t16 * t18 - t17 * t19 + (t4 * t85 + t5 * t84 + t155) * pkin(3), t33, t26, t49, t50, 0, t99 * t86 + (-t103 - t158) * t90 + t141, t99 * t90 + (t103 - t152) * t86 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, t50, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 * t73 * t90, t134 * t73, t86 * t74, t142, qJDD(5), t102 * t86 + t125 * t90, t102 * t90 - t125 * t86;];
tau_reg = t1;
