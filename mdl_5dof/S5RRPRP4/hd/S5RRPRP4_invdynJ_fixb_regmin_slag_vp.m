% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:58
% EndTime: 2019-12-31 19:53:00
% DurationCPUTime: 0.89s
% Computational Cost: add. (1135->209), mult. (1394->220), div. (0->0), fcn. (651->8), ass. (0->130)
t104 = qJD(4) * pkin(4) - qJD(5);
t137 = pkin(1) * qJD(1);
t78 = cos(qJ(2));
t116 = t78 * t137;
t100 = qJD(3) - t116;
t70 = qJD(1) + qJD(2);
t80 = -pkin(2) - pkin(7);
t23 = t80 * t70 + t100;
t77 = cos(qJ(4));
t149 = t77 * t23;
t14 = -t104 - t149;
t127 = qJD(4) * qJ(5);
t74 = sin(qJ(4));
t150 = t74 * t23;
t16 = t127 + t150;
t123 = qJDD(4) * qJ(5);
t113 = qJD(2) * t137;
t133 = pkin(1) * qJDD(1);
t75 = sin(qJ(2));
t141 = -t75 * t113 + t78 * t133;
t112 = qJDD(3) - t141;
t69 = qJDD(1) + qJDD(2);
t15 = t80 * t69 + t112;
t9 = t74 * t15;
t4 = t123 + t9 + (qJD(5) + t149) * qJD(4);
t10 = t77 * t15;
t128 = qJDD(4) * pkin(4);
t130 = qJD(4) * t74;
t5 = t23 * t130 + qJDD(5) - t10 - t128;
t84 = t4 * t74 - t5 * t77 + (t14 * t74 + t16 * t77) * qJD(4);
t165 = t74 * pkin(4) - t77 * qJ(5);
t162 = qJ(3) + t165;
t167 = t162 * t69;
t166 = t162 * t70;
t73 = qJ(1) + qJ(2);
t64 = cos(t73);
t54 = g(2) * t64;
t63 = sin(t73);
t55 = g(1) * t63;
t142 = t55 - t54;
t82 = t142 - t84;
t151 = t64 * t77;
t163 = g(2) * t151 + g(3) * t74;
t101 = -g(1) * t64 - g(2) * t63;
t81 = qJD(4) ^ 2;
t148 = t80 * t81;
t103 = -t77 * qJD(5) + qJD(3);
t129 = qJD(4) * t77;
t22 = pkin(4) * t129 + t74 * t127 + t103;
t161 = -(-t22 + t116) * t70 + t167 - t148;
t68 = t70 ^ 2;
t117 = t75 * t137;
t13 = t117 + t166;
t144 = t78 * t113 + t75 * t133;
t97 = pkin(4) * t77 + qJ(5) * t74;
t2 = t167 + (t97 * qJD(4) + t103) * t70 + t144;
t160 = t13 * t129 + t2 * t74;
t76 = sin(qJ(1));
t159 = g(1) * t76;
t157 = t69 * pkin(2);
t156 = t75 * pkin(1);
t155 = t13 * t70;
t135 = t70 * qJ(3);
t30 = t117 + t135;
t154 = t30 * t70;
t58 = -t78 * pkin(1) - pkin(2);
t47 = -pkin(7) + t58;
t153 = t47 * t81;
t152 = t63 * t77;
t44 = t77 * t69;
t136 = t69 * qJ(3);
t17 = t70 * qJD(3) + t136 + t144;
t147 = t30 * t129 + t17 * t74;
t132 = qJD(2) * t75;
t119 = pkin(1) * t132;
t62 = qJDD(4) * t77;
t146 = t119 * t129 + t47 * t62;
t145 = -g(1) * t151 - g(2) * t152;
t143 = t64 * pkin(2) + t63 * qJ(3);
t140 = t68 + t81;
t71 = t74 ^ 2;
t72 = t77 ^ 2;
t139 = t71 - t72;
t138 = t71 + t72;
t131 = qJD(2) * t78;
t126 = qJDD(4) * t47;
t125 = qJDD(4) * t74;
t124 = qJDD(4) * t80;
t122 = t13 * t130 - t145;
t121 = t10 + t163;
t120 = t17 * t77 + t145;
t118 = pkin(1) * t131;
t115 = t70 * t132;
t114 = t70 * t129;
t49 = t64 * qJ(3);
t111 = -t63 * pkin(2) + t49;
t110 = t138 * t69;
t109 = -t141 - t142;
t108 = t144 + t101;
t106 = pkin(1) * t115;
t105 = t70 * t117;
t99 = t155 + t55;
t98 = -t154 - t55;
t96 = t14 * t77 - t16 * t74;
t20 = t112 - t157;
t92 = -t117 + t135;
t91 = t64 * pkin(7) + t165 * t63 + t143;
t18 = t22 + t118;
t28 = t162 + t156;
t90 = -t18 * t70 - t28 * t69 + t153;
t42 = qJD(3) + t118;
t50 = qJ(3) + t156;
t89 = t42 * t70 + t50 * t69 - t153;
t88 = (-t117 + t166) * qJD(4);
t87 = t105 - t109;
t85 = g(1) * (t165 * t64 + t80 * t63 + t49);
t83 = t100 * t70 + t136 - t148;
t79 = cos(qJ(1));
t66 = t79 * pkin(1);
t43 = t77 * t124;
t36 = t77 * t68 * t74;
t34 = -t81 * t74 + t62;
t33 = -t81 * t77 - t125;
t29 = -t70 * pkin(2) + t100;
t26 = t97 * t70;
t25 = -t140 * t74 + t62;
t24 = t140 * t77 + t125;
t21 = -0.2e1 * t74 * t114 + t72 * t69;
t8 = 0.2e1 * t139 * t70 * qJD(4) - 0.2e1 * t74 * t44;
t1 = [qJDD(1), -g(2) * t79 + t159, g(1) * t79 + g(2) * t76, t69, (t69 * t78 - t115) * pkin(1) - t109, (-t70 * t131 - t69 * t75) * pkin(1) - t108, t106 + qJDD(3) + (-pkin(2) + t58) * t69 + t109, (qJD(3) + t42) * t70 + (qJ(3) + t50) * t69 + t108, t17 * t50 + t30 * t42 + t20 * t58 + t29 * t119 - g(1) * (-t76 * pkin(1) + t111) - g(2) * (t66 + t143), t21, t8, t34, t33, 0, t50 * t114 + (t101 + t89) * t74 + t146 + t147, t89 * t77 + (-t126 + (-t50 * t70 - t119 - t30) * qJD(4)) * t74 + t120, t28 * t114 + (t101 - t90) * t74 + t146 + t160, -t138 * t106 - t47 * t110 + t82, (t126 + (t28 * t70 + t119) * qJD(4)) * t74 + (-t2 + t90) * t77 + t122, t2 * t28 + t13 * t18 - t85 - g(2) * (t66 + t91) + (-t96 * t132 + t159) * pkin(1) + t84 * t47; 0, 0, 0, t69, t87, t70 * t116 - t108, qJDD(3) - t87 - 0.2e1 * t157, 0.2e1 * t136 + (0.2e1 * qJD(3) - t116) * t70 + t108, t17 * qJ(3) + t30 * qJD(3) - t20 * pkin(2) - g(1) * t111 - g(2) * t143 + (-t29 * t75 - t30 * t78) * t137, t21, t8, t34, t33, 0, t43 + t92 * t129 + (t101 + t83) * t74 + t147, t83 * t77 + (-t124 + (-t30 - t92) * qJD(4)) * t74 + t120, t43 + t77 * t88 + (t101 + t161) * t74 + t160, t138 * t105 - t80 * t110 + t82, (t88 + t124) * t74 + (-t161 - t2) * t77 + t122, t2 * t162 + t13 * t22 - t85 - g(2) * t91 + (-t13 * t78 + t96 * t75) * t137 + t84 * t80; 0, 0, 0, 0, 0, 0, t69, -t68, -t142 + t20 - t154, 0, 0, 0, 0, 0, t25, -t24, t25, -t110, t24, -t155 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t139 * t68, t44, -t74 * t69, qJDD(4), t98 * t77 + t121, g(3) * t77 - t9 + (-t98 - t54) * t74, -g(1) * t152 + 0.2e1 * t128 - qJDD(5) + (-t13 * t77 - t26 * t74) * t70 + t121, -t97 * t69 + ((t16 - t127) * t77 + (t104 + t14) * t74) * t70, 0.2e1 * t123 + 0.2e1 * qJD(4) * qJD(5) + t9 + (t26 * t70 - g(3)) * t77 + (-t99 + t54) * t74, t4 * qJ(5) - t5 * pkin(4) - t13 * t26 - t14 * t150 + g(3) * t165 + (qJD(5) - t149) * t16 - t142 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t36, t44, -t72 * t68 - t81, -t16 * qJD(4) + t99 * t77 - t163 + t5;];
tau_reg = t1;
