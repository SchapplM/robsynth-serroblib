% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:56
% EndTime: 2019-03-08 20:21:01
% DurationCPUTime: 1.77s
% Computational Cost: add. (2153->299), mult. (4747->405), div. (0->0), fcn. (3148->8), ass. (0->161)
t85 = sin(qJ(4));
t152 = qJD(2) * t85;
t77 = qJD(5) + t152;
t119 = qJD(5) * t85 + qJD(2);
t88 = cos(qJ(4));
t151 = qJD(2) * t88;
t87 = cos(qJ(5));
t131 = t87 * t151;
t144 = t87 * qJD(4);
t84 = sin(qJ(5));
t146 = qJD(5) * t84;
t130 = t88 * t146;
t99 = t85 * t144 + t130;
t42 = t99 * qJD(2) - qJD(5) * t144;
t165 = t88 * t42;
t176 = t77 * t87;
t177 = t77 * t84;
t150 = qJD(4) * t84;
t65 = t131 + t150;
t201 = ((-t65 + t131) * t85 + t88 * t176) * qJD(4) - t119 * t177 - t165;
t142 = qJD(2) * qJD(4);
t125 = t88 * t142;
t118 = pkin(5) * t125;
t145 = qJD(5) * t87;
t82 = sin(pkin(6));
t153 = qJD(2) * t82;
t86 = sin(qJ(2));
t133 = t86 * t153;
t148 = qJD(4) * t88;
t154 = qJD(1) * t85;
t155 = qJD(1) * t82;
t89 = cos(qJ(2));
t127 = t89 * t155;
t109 = qJD(3) - t127;
t90 = -pkin(2) - pkin(8);
t56 = t90 * qJD(2) + t109;
t83 = cos(pkin(6));
t21 = t56 * t148 + (-qJD(4) * t83 + t133) * t154;
t173 = t83 * t88;
t76 = qJD(1) * t173;
t38 = t85 * t56 + t76;
t29 = qJD(4) * pkin(9) + t38;
t116 = pkin(4) * t88 + pkin(9) * t85;
t60 = t116 * qJD(4) + qJD(3);
t39 = (t60 + t127) * qJD(2);
t128 = t86 * t155;
t69 = pkin(4) * t85 - pkin(9) * t88 + qJ(3);
t47 = t69 * qJD(2) + t128;
t124 = -t29 * t145 - t47 * t146 - t84 * t21 + t87 * t39;
t2 = -t118 - t124;
t10 = t29 * t87 + t47 * t84;
t7 = qJ(6) * t77 + t10;
t200 = -t7 * t77 + t2;
t169 = t85 * t90;
t160 = t87 * t169 + t84 * t69;
t171 = t84 * t86;
t199 = -t87 * t60 - (t85 * t171 - t87 * t89) * t155;
t168 = t86 * t87;
t198 = (t85 * t168 + t84 * t89) * t155 - t88 * t90 * t144 - t69 * t145 - t84 * t60;
t147 = qJD(5) * t65;
t172 = t84 * t85;
t73 = t142 * t172;
t43 = -t73 + t147;
t37 = -t83 * t154 + t56 * t88;
t195 = t65 ^ 2;
t122 = t88 * t133;
t149 = qJD(4) * t85;
t161 = -qJD(4) * t76 - t56 * t149;
t22 = -qJD(1) * t122 - t161;
t3 = pkin(5) * t43 + qJ(6) * t42 - qJD(6) * t65 + t22;
t194 = t3 * t84;
t193 = t3 * t87;
t191 = -qJ(6) * t148 - (-t90 * t146 + qJD(6)) * t85 + t198;
t28 = -qJD(4) * pkin(4) - t37;
t63 = t84 * t151 - t144;
t11 = pkin(5) * t63 - qJ(6) * t65 + t28;
t190 = t11 * t65;
t189 = t22 * t84;
t188 = t22 * t87;
t187 = t28 * t84;
t186 = t28 * t87;
t185 = t42 * t84;
t184 = t42 * t85;
t183 = t43 * t85;
t181 = t63 * t77;
t180 = t65 * t63;
t179 = t65 * t77;
t143 = qJD(2) * qJ(3);
t68 = t128 + t143;
t178 = t68 * t89;
t175 = t82 * t89;
t92 = qJD(2) ^ 2;
t174 = t82 * t92;
t170 = t84 * t90;
t166 = t87 * t69;
t126 = -pkin(5) + t170;
t164 = t160 * qJD(5) + t126 * t148 + t199;
t110 = pkin(5) * t84 - qJ(6) * t87;
t163 = t84 * qJD(6) - t77 * t110 + t38;
t67 = t116 * qJD(2);
t162 = t87 * t37 + t84 * t67;
t81 = t88 ^ 2;
t159 = t85 ^ 2 - t81;
t91 = qJD(4) ^ 2;
t158 = -t91 - t92;
t157 = qJD(2) * pkin(2);
t9 = -t29 * t84 + t47 * t87;
t156 = qJD(6) - t9;
t141 = pkin(9) * t177;
t140 = pkin(9) * t176;
t138 = t86 * t174;
t137 = t89 * t174;
t136 = -t47 * t145 - t87 * t21 - t84 * t39;
t134 = pkin(9) * t148;
t132 = t89 * t153;
t129 = t77 * t145;
t123 = t88 * t128;
t121 = t65 * t128;
t120 = t85 * t133;
t117 = qJ(6) * t125;
t115 = -t68 + t128;
t6 = -pkin(5) * t77 + t156;
t113 = t6 * t87 - t7 * t84;
t112 = t6 * t84 + t7 * t87;
t111 = pkin(5) * t87 + qJ(6) * t84;
t108 = -t37 * t84 + t67 * t87;
t107 = qJD(2) * t81 - t77 * t85;
t55 = -t85 * t175 + t173;
t106 = t82 * t168 - t55 * t84;
t34 = t82 * t171 + t55 * t87;
t54 = t88 * t175 + t83 * t85;
t105 = t110 - t90;
t104 = -t11 * t85 + t134;
t103 = t28 * t85 - t134;
t102 = t10 * t77 + t124;
t101 = t115 * qJD(2);
t100 = t29 * t146 + t136;
t98 = t115 - t143;
t32 = t55 * qJD(4) - t122;
t31 = -t54 * qJD(4) + t120;
t4 = t34 * qJD(5) - t87 * t132 + t31 * t84;
t97 = t106 * t125 + t32 * t63 - t4 * t77 + t54 * t43;
t61 = (qJD(3) + t127) * qJD(2);
t96 = t109 * qJD(2) - t90 * t91 + t61;
t1 = qJD(6) * t77 - t100 + t117;
t95 = t113 * qJD(5) + t1 * t87 + t2 * t84;
t5 = t106 * qJD(5) + t84 * t132 + t31 * t87;
t94 = t34 * t125 - t32 * t65 + t42 * t54 + t5 * t77;
t93 = t63 * t149 + (-t43 - t73) * t88 + (-t119 * t87 - t84 * t148) * t77;
t70 = -pkin(4) - t111;
t62 = t109 - t157;
t48 = t105 * t88;
t46 = t63 * t123;
t41 = t126 * t85 - t166;
t40 = qJ(6) * t85 + t160;
t30 = pkin(5) * t65 + qJ(6) * t63;
t18 = t181 - t42;
t15 = (t111 * qJD(5) - qJD(6) * t87) * t88 - t105 * t149;
t14 = -pkin(5) * t151 - t108;
t13 = qJ(6) * t151 + t162;
t8 = [0, 0, -t138, -t137, t138, t137 (t61 * t86 + (t178 + (t62 - t127) * t86) * qJD(2)) * t82, 0, 0, 0, 0, 0, t85 * t137 + (-t32 + t122) * qJD(4), t88 * t137 + (-t31 - t120) * qJD(4), 0, 0, 0, 0, 0, t97, -t94, t97, t106 * t42 - t34 * t43 + t4 * t65 - t5 * t63, t94, t1 * t34 - t106 * t2 + t11 * t32 + t3 * t54 + t4 * t6 + t5 * t7; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t61 * qJ(3) + t68 * qJD(3) + (-t178 + (-t62 - t157) * t86) * t155, -0.2e1 * t85 * t125, 0.2e1 * t159 * t142, -t91 * t85, -t91 * t88, 0, -t98 * t148 + t96 * t85, t98 * t149 + t96 * t88, -t87 * t165 - t99 * t65 (t63 * t87 + t65 * t84) * t149 + (t185 - t43 * t87 + (t63 * t84 - t65 * t87) * qJD(5)) * t88, -t77 * t130 - t184 + (t107 * t87 + t65 * t88) * qJD(4), -t88 * t129 - t183 + (-t107 * t84 - t63 * t88) * qJD(4) (t77 + t152) * t148, t46 + (-t146 * t69 - t199) * t77 + (-t90 * t129 + (t63 * t90 - t187) * qJD(4) + t124) * t85 + (t28 * t145 + t189 - t90 * t43 + (-t77 * t170 + (-t169 * t84 + t166) * qJD(2) + t9) * qJD(4)) * t88, t198 * t77 + ((t77 * t90 + t29) * t146 + (t65 * t90 - t186) * qJD(4) + t136) * t85 + (t121 - t28 * t146 + t188 + t90 * t42 + (-qJD(2) * t160 - t10) * qJD(4)) * t88, t15 * t63 + t43 * t48 + t46 + (-t11 * t150 - t2) * t85 - t164 * t77 + (t11 * t145 + t194 + (-qJD(2) * t41 - t6) * qJD(4)) * t88, -t40 * t43 - t41 * t42 + t164 * t65 + t191 * t63 - t113 * t149 + (-qJD(5) * t112 - t1 * t84 + t2 * t87) * t88, -t15 * t65 + t42 * t48 + (t11 * t144 + t1) * t85 - t191 * t77 + (-t121 + t11 * t146 - t193 + (qJD(2) * t40 + t7) * qJD(4)) * t88, t1 * t40 + t2 * t41 + t3 * t48 - t191 * t7 + t164 * t6 + (t15 + t123) * t11; 0, 0, 0, 0, 0, -t92, t101, 0, 0, 0, 0, 0, t158 * t85, t158 * t88, 0, 0, 0, 0, 0, t93, -t201, t93 (t119 * t65 - t148 * t63 - t183) * t87 + (t119 * t63 + t148 * t65 - t184) * t84, t201, t113 * qJD(2) + (qJD(4) * t112 - t3) * t88 + (qJD(4) * t11 + t95) * t85; 0, 0, 0, 0, 0, 0, 0, t88 * t92 * t85, -t159 * t92, 0, 0, 0, t38 * qJD(4) + t88 * t101 + t161, -t115 * t152, t65 * t176 - t185 (-t42 - t181) * t87 + (-t43 - t179) * t84, t129 + (t85 * t176 + (-t65 + t150) * t88) * qJD(2), -t77 * t146 + (-t77 * t172 + (t63 + t144) * t88) * qJD(2), -t77 * t151, -pkin(4) * t43 - t188 - t108 * t77 - t38 * t63 + (-t140 + t187) * qJD(5) + (t103 * t84 - t9 * t88) * qJD(2), pkin(4) * t42 + t189 + t162 * t77 - t38 * t65 + (t141 + t186) * qJD(5) + (t10 * t88 + t103 * t87) * qJD(2), t14 * t77 - t193 + t43 * t70 - t163 * t63 + (t11 * t84 - t140) * qJD(5) + (-t104 * t84 + t6 * t88) * qJD(2), t13 * t63 - t14 * t65 + (t1 + t77 * t6 + (-t43 + t147) * pkin(9)) * t87 + ((qJD(5) * t63 - t42) * pkin(9) + t200) * t84, -t13 * t77 - t194 + t42 * t70 + t163 * t65 + (-t11 * t87 - t141) * qJD(5) + (t104 * t87 - t7 * t88) * qJD(2), pkin(9) * t95 - t11 * t163 - t7 * t13 - t6 * t14 + t3 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, -t63 ^ 2 + t195, t18, t179 - t43, t125, -t28 * t65 + t102, t28 * t63 + t77 * t9 + t100, -t30 * t63 + t102 + 0.2e1 * t118 - t190, pkin(5) * t42 - t43 * qJ(6) + (-t10 + t7) * t65 + (t6 - t156) * t63, 0.2e1 * t117 - t11 * t63 + t30 * t65 + (0.2e1 * qJD(6) - t9) * t77 - t100, -t2 * pkin(5) + t1 * qJ(6) - t6 * t10 - t11 * t30 + t156 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 + t180, t18, -t77 ^ 2 - t195, t190 + t200;];
tauc_reg  = t8;
