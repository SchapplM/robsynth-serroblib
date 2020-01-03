% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:53:01
% EndTime: 2019-12-31 19:53:03
% DurationCPUTime: 1.19s
% Computational Cost: add. (1367->236), mult. (1724->239), div. (0->0), fcn. (790->8), ass. (0->147)
t121 = qJD(4) * pkin(4) - qJD(5);
t156 = pkin(1) * qJD(1);
t85 = cos(qJ(2));
t133 = t85 * t156;
t112 = qJD(3) - t133;
t77 = qJD(1) + qJD(2);
t87 = -pkin(2) - pkin(7);
t26 = t87 * t77 + t112;
t84 = cos(qJ(4));
t168 = t84 * t26;
t16 = -t121 - t168;
t146 = qJD(4) * qJ(5);
t81 = sin(qJ(4));
t170 = t81 * t26;
t18 = t146 + t170;
t130 = qJD(2) * t156;
t152 = pkin(1) * qJDD(1);
t82 = sin(qJ(2));
t160 = -t82 * t130 + t85 * t152;
t129 = qJDD(3) - t160;
t76 = qJDD(1) + qJDD(2);
t17 = t87 * t76 + t129;
t10 = t81 * t17;
t142 = qJDD(4) * qJ(5);
t4 = t142 + t10 + (qJD(5) + t168) * qJD(4);
t11 = t84 * t17;
t147 = qJDD(4) * pkin(4);
t149 = qJD(4) * t81;
t5 = t26 * t149 + qJDD(5) - t11 - t147;
t92 = t4 * t81 - t5 * t84 + (t16 * t81 + t18 * t84) * qJD(4);
t78 = t81 ^ 2;
t79 = t84 ^ 2;
t157 = t78 + t79;
t189 = t81 * pkin(4) - t84 * qJ(5);
t185 = qJ(3) + t189;
t191 = t185 * t76;
t190 = t185 * t77;
t80 = qJ(1) + qJ(2);
t71 = cos(t80);
t61 = g(2) * t71;
t70 = sin(t80);
t62 = g(1) * t70;
t161 = t61 - t62;
t90 = -t161 - t92;
t171 = t71 * t84;
t187 = g(2) * t171 + g(3) * t81;
t114 = -g(1) * t71 - g(2) * t70;
t127 = t157 * t17;
t186 = t157 * t26;
t134 = t82 * t156;
t122 = t77 * t134;
t31 = t157 * t76;
t184 = t157 * t122 - t87 * t31;
t88 = qJD(4) ^ 2;
t167 = t87 * t88;
t119 = -t84 * qJD(5) + qJD(3);
t148 = qJD(4) * t84;
t25 = pkin(4) * t148 + t81 * t146 + t119;
t183 = -(-t25 + t133) * t77 + t191 - t167;
t75 = t77 ^ 2;
t158 = t78 - t79;
t51 = t84 * t76;
t182 = 0.2e1 * t158 * t77 * qJD(4) - 0.2e1 * t81 * t51;
t14 = t134 + t190;
t109 = pkin(4) * t84 + qJ(5) * t81;
t163 = t85 * t130 + t82 * t152;
t2 = t191 + (t109 * qJD(4) + t119) * t77 + t163;
t181 = t14 * t148 + t2 * t81;
t83 = sin(qJ(1));
t180 = g(1) * t83;
t178 = t76 * pkin(2);
t177 = t82 * pkin(1);
t176 = t14 * t77;
t154 = t77 * qJ(3);
t35 = t134 + t154;
t175 = t35 * t77;
t174 = t35 * t85;
t65 = -t85 * pkin(1) - pkin(2);
t54 = -pkin(7) + t65;
t173 = t54 * t88;
t172 = t70 * t84;
t169 = t81 * t76;
t155 = t76 * qJ(3);
t19 = t77 * qJD(3) + t155 + t163;
t166 = t35 * t148 + t19 * t81;
t151 = qJD(2) * t82;
t136 = pkin(1) * t151;
t69 = qJDD(4) * t84;
t165 = t136 * t148 + t54 * t69;
t164 = -g(1) * t171 - g(2) * t172;
t162 = t71 * pkin(2) + t70 * qJ(3);
t159 = t75 + t88;
t150 = qJD(2) * t85;
t145 = qJDD(4) * t54;
t144 = qJDD(4) * t81;
t143 = qJDD(4) * t87;
t141 = t14 * t149 - t164;
t140 = t11 + t187;
t139 = t19 * t84 + t164;
t58 = t71 * pkin(7);
t138 = t58 + t162;
t86 = cos(qJ(1));
t73 = t86 * pkin(1);
t137 = t73 + t162;
t135 = pkin(1) * t150;
t132 = t77 * t151;
t131 = t77 * t148;
t56 = t71 * qJ(3);
t128 = -t70 * pkin(2) + t56;
t126 = -t160 + t161;
t125 = t163 + t114;
t123 = pkin(1) * t132;
t120 = t81 * t131;
t117 = t161 - t175;
t115 = t87 * t70 + t56;
t113 = -g(2) * t86 + t180;
t111 = t176 + t62;
t110 = -t175 - t62;
t108 = t16 * t84 - t18 * t81;
t49 = qJD(3) + t135;
t57 = qJ(3) + t177;
t106 = t19 * t57 + t35 * t49;
t104 = t19 * qJ(3) + t35 * qJD(3);
t22 = t129 - t178;
t103 = g(1) * t115;
t101 = -t134 + t154;
t100 = t189 * t70 + t138;
t20 = t25 + t135;
t32 = t185 + t177;
t99 = -t20 * t77 - t32 * t76 + t173;
t98 = t49 * t77 + t57 * t76 - t173;
t97 = (-t134 + t190) * qJD(4);
t95 = t122 - t126;
t93 = g(1) * (t189 * t71 + t115);
t91 = t112 * t77 + t155 - t167;
t50 = t84 * t143;
t43 = t84 * t75 * t81;
t41 = -t88 * t81 + t69;
t40 = t88 * t84 + t144;
t34 = -t77 * pkin(2) + t112;
t33 = t158 * t75;
t29 = t109 * t77;
t28 = -t159 * t81 + t69;
t27 = t159 * t84 + t144;
t24 = t79 * t76 - 0.2e1 * t120;
t23 = t78 * t76 + 0.2e1 * t120;
t1 = [0, 0, 0, 0, 0, qJDD(1), t113, g(1) * t86 + g(2) * t83, 0, 0, 0, 0, 0, 0, 0, t76, (t76 * t85 - t132) * pkin(1) - t126, (-t77 * t150 - t76 * t82) * pkin(1) - t125, 0, (t113 + (t82 ^ 2 + t85 ^ 2) * t152) * pkin(1), t76, 0, 0, 0, 0, 0, 0, t123 + qJDD(3) + (-pkin(2) + t65) * t76 + t126, (qJD(3) + t49) * t77 + (qJ(3) + t57) * t76 + t125, t22 * t65 + t34 * t136 - g(1) * (-t83 * pkin(1) + t128) - g(2) * t137 + t106, t24, t182, t41, t23, -t40, 0, t57 * t131 + (t114 + t98) * t81 + t165 + t166, t98 * t84 + (-t145 + (-t57 * t77 - t136 - t35) * qJD(4)) * t81 + t139, -t161 + t157 * (-t54 * t76 - t123 - t17), -t103 - g(2) * (t58 + t137) + t54 * t127 + (t151 * t186 + t180) * pkin(1) + t106, t24, t41, -t182, 0, t40, t23, t32 * t131 + (t114 - t99) * t81 + t165 + t181, -t123 * t157 - t31 * t54 + t90, (t145 + (t32 * t77 + t136) * qJD(4)) * t81 + (-t2 + t99) * t84 + t141, t2 * t32 + t14 * t20 - t93 - g(2) * (t100 + t73) + (-t108 * t151 + t180) * pkin(1) + t92 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t95, t77 * t133 - t125, 0, 0, t76, 0, 0, 0, 0, 0, 0, qJDD(3) - t95 - 0.2e1 * t178, 0.2e1 * t155 + (0.2e1 * qJD(3) - t133) * t77 + t125, -t22 * pkin(2) - g(1) * t128 - g(2) * t162 + (-t34 * t82 - t174) * t156 + t104, t24, t182, t41, t23, -t40, 0, t50 + t101 * t148 + (t114 + t91) * t81 + t166, t91 * t84 + (-t143 + (-t101 - t35) * qJD(4)) * t81 + t139, -t161 + t184 - t127, -t103 - g(2) * t138 + t87 * t127 + (-t82 * t186 - t174) * t156 + t104, t24, t41, -t182, 0, t40, t23, t50 + t84 * t97 + (t114 + t183) * t81 + t181, t90 + t184, (t97 + t143) * t81 + (-t183 - t2) * t84 + t141, t2 * t185 + t14 * t25 - t93 - g(2) * t100 + (t108 * t82 - t14 * t85) * t156 + t92 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t75, t117 + t22, 0, 0, 0, 0, 0, 0, t28, -t27, -t31, t127 + t117, 0, 0, 0, 0, 0, 0, t28, -t31, t27, -t176 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t33, t51, -t43, -t169, qJDD(4), t110 * t84 + t140, g(3) * t84 - t10 + (-t110 - t61) * t81, 0, 0, t43, t51, t33, qJDD(4), t169, -t43, -g(1) * t172 + 0.2e1 * t147 - qJDD(5) + (-t14 * t84 - t29 * t81) * t77 + t140, -t109 * t76 + ((t18 - t146) * t84 + (t121 + t16) * t81) * t77, 0.2e1 * t142 + 0.2e1 * qJD(4) * qJD(5) + t10 + (t29 * t77 - g(3)) * t84 + (-t111 + t61) * t81, t4 * qJ(5) - t5 * pkin(4) - t14 * t29 - t16 * t170 + g(3) * t189 + (qJD(5) - t168) * t18 + t161 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t43, t51, -t79 * t75 - t88, -t18 * qJD(4) + t111 * t84 - t187 + t5;];
tau_reg = t1;
