% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:44
% EndTime: 2019-12-31 19:49:47
% DurationCPUTime: 0.86s
% Computational Cost: add. (1365->208), mult. (2182->244), div. (0->0), fcn. (1233->12), ass. (0->140)
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t115 = t98 * pkin(4) + t95 * qJ(5);
t111 = -pkin(3) - t115;
t153 = pkin(1) * qJD(2);
t131 = qJD(1) * t153;
t96 = sin(qJ(2));
t142 = qJDD(1) * t96;
t99 = cos(qJ(2));
t183 = pkin(1) * t142 + t99 * t131;
t169 = t99 * pkin(1);
t79 = qJDD(1) * t169;
t88 = qJDD(1) + qJDD(2);
t33 = t88 * pkin(2) - t96 * t131 + t79;
t93 = sin(pkin(8));
t94 = cos(pkin(8));
t16 = t183 * t94 + t93 * t33;
t185 = t88 * pkin(7) + qJD(3) * qJD(4) + t16;
t154 = pkin(1) * qJD(1);
t133 = t99 * t154;
t89 = qJD(1) + qJD(2);
t45 = t89 * pkin(2) + t133;
t134 = t96 * t154;
t61 = t94 * t134;
t30 = t93 * t45 + t61;
t26 = t89 * pkin(7) + t30;
t162 = t95 * t26;
t18 = t98 * qJD(3) - t162;
t184 = qJD(5) - t18;
t137 = t95 * qJDD(3) + t185 * t98;
t138 = qJDD(4) * qJ(5);
t3 = t138 + (qJD(5) - t162) * qJD(4) + t137;
t147 = qJD(4) * t98;
t124 = -t98 * qJDD(3) + t26 * t147 + t185 * t95;
t146 = qJDD(4) * pkin(4);
t178 = qJDD(5) - t146;
t4 = t124 + t178;
t182 = t3 * t98 + t4 * t95;
t90 = t95 ^ 2;
t91 = t98 ^ 2;
t155 = t90 + t91;
t181 = t155 * t89;
t14 = -qJD(4) * pkin(4) + t184;
t143 = qJD(4) * qJ(5);
t19 = t95 * qJD(3) + t98 * t26;
t17 = t19 + t143;
t92 = qJ(1) + qJ(2);
t84 = sin(t92);
t85 = cos(t92);
t180 = g(1) * t84 - g(2) * t85;
t83 = pkin(8) + t92;
t68 = sin(t83);
t69 = cos(t83);
t179 = g(1) * t69 + g(2) * t68;
t177 = pkin(2) * t84;
t176 = g(1) * t68;
t173 = g(2) * t69;
t171 = t94 * pkin(2);
t37 = t93 * t133 + t61;
t168 = t37 * t89;
t167 = t68 * t95;
t166 = t69 * t95;
t165 = t89 * t95;
t164 = t93 * t96;
t163 = t94 * t96;
t161 = t98 * t88;
t144 = t95 * qJD(5);
t148 = qJD(4) * t95;
t41 = pkin(4) * t148 - t98 * t143 - t144;
t160 = t41 - t37;
t159 = g(1) * t167 - g(2) * t166;
t78 = pkin(2) + t169;
t158 = pkin(1) * t163 + t93 * t78;
t157 = g(1) * t85 + g(2) * t84;
t156 = t90 - t91;
t101 = qJD(4) ^ 2;
t36 = pkin(7) + t158;
t152 = t101 * t36;
t70 = t93 * pkin(2) + pkin(7);
t151 = t101 * t70;
t149 = qJD(4) * t89;
t145 = t17 * qJD(4);
t141 = qJDD(4) * t36;
t140 = qJDD(4) * t70;
t87 = t89 ^ 2;
t136 = t95 * t87 * t98;
t60 = t93 * t134;
t39 = t94 * t133 - t60;
t53 = t98 * t176;
t135 = t39 * t148 + t98 * t168 + t53;
t15 = -t183 * t93 + t94 * t33;
t10 = -t88 * pkin(3) - t15;
t130 = -t10 - t173;
t129 = t155 * t88;
t128 = t18 + t162;
t126 = -pkin(1) * t164 + t94 * t78;
t27 = t111 - t126;
t40 = (t94 * t99 - t164) * t153;
t127 = t27 * t89 - t40;
t29 = t94 * t45 - t60;
t25 = -t89 * pkin(3) - t29;
t125 = t10 * t95 + t25 * t147 - t159;
t122 = qJD(1) * (-qJD(2) + t89);
t121 = qJD(2) * (-qJD(1) - t89);
t119 = t79 + t180;
t77 = pkin(2) * t85;
t118 = t68 * pkin(7) - t111 * t69 + t77;
t97 = sin(qJ(1));
t117 = -t97 * pkin(1) - t177;
t114 = pkin(4) * t95 - qJ(5) * t98;
t113 = t14 * t95 + t17 * t98;
t71 = -pkin(3) - t171;
t112 = t71 * t88 + t151;
t44 = t111 - t171;
t5 = (t114 * qJD(4) - t144) * t89 + t111 * t88 - t15;
t110 = -t44 * t88 - t151 - t5;
t38 = (t93 * t99 + t163) * t153;
t35 = -pkin(3) - t126;
t109 = t35 * t88 + t38 * t89 + t152;
t108 = g(1) * t166 + g(2) * t167 - g(3) * t98 - t124;
t107 = t111 * t176;
t20 = t38 + t41;
t106 = -t20 * t89 - t27 * t88 - t152 - t5;
t105 = -t141 + (t35 * t89 - t40) * qJD(4);
t104 = t19 * qJD(4) + t108;
t103 = t14 * t147 - t95 * t145 - t179 + t182;
t102 = (t14 * t98 - t17 * t95) * qJD(4) + t182;
t100 = cos(qJ(1));
t86 = t100 * pkin(1);
t67 = t95 * t88;
t63 = t69 * pkin(7);
t48 = qJDD(4) * t98 - t101 * t95;
t47 = qJDD(4) * t95 + t101 * t98;
t42 = t114 * t89;
t34 = 0.2e1 * t147 * t165 + t90 * t88;
t24 = -0.2e1 * t156 * t149 + 0.2e1 * t95 * t161;
t21 = t25 * t148;
t13 = t111 * t89 - t29;
t9 = t13 * t148;
t1 = [qJDD(1), g(1) * t97 - g(2) * t100, g(1) * t100 + g(2) * t97, t88, (t96 * t121 + t88 * t99) * pkin(1) + t119, ((-qJDD(1) - t88) * t96 + t99 * t121) * pkin(1) + t157, t16 * t158 + t30 * t40 + t15 * t126 - t29 * t38 - g(1) * t117 - g(2) * (t77 + t86), t34, t24, t47, t48, 0, t21 + t53 + t105 * t95 + (-t109 + t130) * t98, t105 * t98 + t109 * t95 + t125, t53 + t9 + (t127 * qJD(4) - t141) * t95 + (t106 - t173) * t98, t36 * t129 + t40 * t181 + t103, (t141 + (-t127 - t13) * qJD(4)) * t98 + t106 * t95 + t159, t5 * t27 + t13 * t20 - g(1) * (t117 + t63) - g(2) * (t86 + t118) - t107 + t113 * t40 + t102 * t36; 0, 0, 0, t88, t96 * pkin(1) * t122 + t119, (t99 * t122 - t142) * pkin(1) + t157, t29 * t37 - t30 * t39 + (t15 * t94 + t16 * t93 + t180) * pkin(2), t34, t24, t47, t48, 0, t21 + (t71 * t149 - t140) * t95 + (-t112 + t130) * t98 + t135, (-t140 + (t71 * t89 + t39) * qJD(4)) * t98 + (t112 - t168) * t95 + t125, t9 + (t44 * t149 - t140) * t95 + (-t41 * t89 + t110 - t173) * t98 + t135, t70 * t129 - t39 * t181 + t103, (t140 + (-t44 * t89 - t13 - t39) * qJD(4)) * t98 + (-t160 * t89 + t110) * t95 + t159, t5 * t44 - g(1) * (t63 - t177) - g(2) * t118 - t107 - t113 * t39 + t160 * t13 + t102 * t70; 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, t48, -t47, t48, 0, t47, t113 * qJD(4) + t3 * t95 - t4 * t98 - g(3); 0, 0, 0, 0, 0, 0, 0, -t136, t156 * t87, t67, t161, qJDD(4), -t25 * t165 + t104, g(3) * t95 + t128 * qJD(4) + (-t25 * t89 + t179) * t98 - t137, 0.2e1 * t146 - qJDD(5) + (-t13 * t95 + t42 * t98) * t89 + t104, -t114 * t88, 0.2e1 * t138 + (t42 * t89 - g(3)) * t95 + (t13 * t89 - t179) * t98 + (0.2e1 * qJD(5) - t128) * qJD(4) + t137, -t4 * pkin(4) - g(3) * t115 + t3 * qJ(5) + t179 * t114 - t13 * t42 - t14 * t19 + t184 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t136, t67, -t90 * t87 - t101, t13 * t165 - t108 - t145 + t178;];
tau_reg = t1;
