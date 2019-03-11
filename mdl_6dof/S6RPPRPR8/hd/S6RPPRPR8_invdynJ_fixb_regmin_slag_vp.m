% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:14
% EndTime: 2019-03-09 01:56:19
% DurationCPUTime: 1.96s
% Computational Cost: add. (2352->324), mult. (4646->377), div. (0->0), fcn. (3318->10), ass. (0->174)
t114 = cos(pkin(9));
t121 = cos(qJ(4));
t189 = t121 * t114;
t173 = qJD(1) * t189;
t118 = sin(qJ(4));
t113 = sin(pkin(9));
t183 = qJD(1) * t113;
t89 = t118 * t183;
t72 = -t89 + t173;
t64 = qJD(6) + t72;
t79 = t113 * t121 + t114 * t118;
t133 = -qJD(4) * t89 + qJDD(1) * t79;
t233 = t133 + qJD(4) * (t72 + t173);
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t228 = g(1) * t119 - g(2) * t122;
t109 = pkin(9) + qJ(4);
t98 = sin(t109);
t99 = cos(t109);
t131 = -g(3) * t99 - t228 * t98;
t181 = qJD(4) * t121;
t182 = qJD(4) * t118;
t116 = -pkin(1) - qJ(3);
t221 = -qJD(1) * qJD(3) + qJDD(1) * t116;
t80 = qJDD(2) + t221;
t165 = -pkin(7) * qJDD(1) + t80;
t52 = t165 * t113;
t53 = t165 * t114;
t87 = t116 * qJD(1) + qJD(2);
t170 = -pkin(7) * qJD(1) + t87;
t60 = t170 * t113;
t61 = t170 * t114;
t166 = -t118 * t53 - t121 * t52 - t61 * t181 + t60 * t182;
t232 = t131 - t166;
t175 = t114 * qJDD(1);
t176 = t113 * qJDD(1);
t149 = t118 * t176 - t121 * t175;
t70 = t79 * qJD(1);
t195 = qJD(4) * t70;
t38 = t149 + t195;
t204 = -t118 * t60 + t121 * t61;
t231 = qJD(5) - t204;
t230 = t64 ^ 2;
t184 = t113 ^ 2 + t114 ^ 2;
t229 = t184 * t87;
t155 = g(1) * t122 + g(2) * t119;
t110 = qJDD(1) * qJ(2);
t111 = qJD(1) * qJD(2);
t227 = t110 + t111;
t84 = qJDD(3) + t227;
t134 = -t155 + t84;
t74 = -t113 * t181 - t114 * t182;
t78 = t113 * t118 - t189;
t203 = t74 * qJD(4) - t78 * qJDD(4);
t226 = qJD(1) * t70 - t203;
t187 = pkin(5) * t72 + t231;
t215 = pkin(5) * t70;
t31 = t118 * t61 + t121 * t60;
t26 = -qJD(4) * qJ(5) - t31;
t16 = -t26 - t215;
t217 = pkin(4) + pkin(8);
t37 = -qJDD(6) + t38;
t224 = t217 * t37 + (t16 - t31 + t215) * t64;
t205 = -pkin(7) + t116;
t81 = t205 * t113;
t82 = t205 * t114;
t41 = t118 * t82 + t121 * t81;
t25 = -qJD(3) * t78 + qJD(4) * t41;
t40 = t118 * t81 - t121 * t82;
t223 = -qJD(4) * t25 - qJDD(4) * t40 - t155 * t98;
t24 = qJD(3) * t79 - t82 * t181 + t182 * t81;
t222 = qJD(4) * t24 - qJDD(4) * t41 - t155 * t99;
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t172 = t114 * t181;
t39 = qJD(1) * t172 + t133;
t47 = qJD(4) * t120 + t117 * t70;
t12 = qJD(6) * t47 + t117 * qJDD(4) - t120 * t39;
t220 = t70 ^ 2;
t219 = t72 ^ 2;
t218 = 0.2e1 * t111;
t216 = pkin(4) * t39;
t214 = g(3) * t98;
t101 = t113 * pkin(3);
t212 = t16 * t79;
t93 = qJ(2) + t101;
t158 = qJ(5) * t78 + t93;
t22 = t217 * t79 + t158;
t211 = t22 * t37;
t45 = qJD(4) * t117 - t120 * t70;
t210 = t45 * t64;
t209 = t45 * t70;
t208 = t47 * t64;
t207 = t47 * t70;
t206 = t72 * t70;
t202 = qJ(5) * t38;
t179 = qJD(6) * t120;
t180 = qJD(6) * t117;
t11 = -qJD(4) * t180 + t120 * qJDD(4) + t117 * t39 + t70 * t179;
t201 = t11 * t120;
t200 = t117 * t37;
t32 = t120 * t37;
t199 = t70 * qJ(5);
t198 = pkin(1) * qJDD(1);
t196 = qJD(4) * t31;
t194 = qJDD(4) * pkin(4);
t193 = t117 * t119;
t192 = t117 * t122;
t191 = t119 * t120;
t190 = t120 * t122;
t186 = t122 * pkin(1) + t119 * qJ(2);
t97 = qJD(1) * qJ(2) + qJD(3);
t177 = qJDD(4) * qJ(5);
t174 = t79 * t179;
t83 = pkin(3) * t183 + t97;
t171 = g(2) * t186;
t13 = -t217 * qJD(4) + t187;
t95 = pkin(3) * t176;
t77 = t95 + t84;
t130 = -qJD(5) * t72 + t202 + t77;
t4 = t217 * t39 + t130;
t169 = qJD(6) * t13 + t4;
t143 = -qJ(5) * t72 + t83;
t17 = t217 * t70 + t143;
t147 = t118 * t52 - t121 * t53 + t60 * t181 + t61 * t182;
t140 = qJDD(5) + t147;
t2 = -pkin(5) * t38 - t217 * qJDD(4) + t140;
t168 = -qJD(6) * t17 + t2;
t167 = t184 * t80;
t164 = t117 * t64;
t163 = -qJD(6) * t78 + qJD(1);
t161 = qJDD(2) - t198;
t156 = -pkin(4) * t98 + qJ(5) * t99;
t75 = -t113 * t182 + t172;
t152 = t11 * t79 + t47 * t75;
t151 = -t37 * t79 + t64 * t75;
t150 = t38 * t78 + t72 * t74;
t6 = t117 * t13 + t120 * t17;
t144 = -qJD(4) * t75 - qJDD(4) * t79;
t142 = -t164 * t64 - t32;
t141 = -qJ(5) * t74 + qJD(5) * t78 + qJD(2);
t28 = -pkin(5) * t78 + t40;
t7 = -qJD(4) * qJD(5) + t166 - t177;
t3 = -pkin(5) * t39 - t7;
t138 = t16 * t75 + t28 * t37 + t3 * t79;
t137 = -t156 + t101;
t136 = -t230 * t120 + t200;
t135 = qJD(1) * t72 - t144;
t23 = -qJD(4) * pkin(4) + t231;
t9 = t140 - t194;
t129 = -t23 * t74 - t26 * t75 - t7 * t79 + t78 * t9 - t228;
t128 = t227 + t134;
t127 = -t228 * t99 - t147 + t214;
t27 = pkin(4) * t70 + t143;
t126 = t27 * t72 + qJDD(5) - t127;
t125 = t3 + (t217 * t64 + t199) * t64 + t131;
t124 = qJD(1) ^ 2;
t115 = -pkin(7) - qJ(3);
t103 = t122 * qJ(2);
t68 = -t99 * t193 + t190;
t67 = -t99 * t191 - t192;
t66 = -t99 * t192 - t191;
t65 = -t99 * t190 + t193;
t36 = pkin(4) * t72 + t199;
t35 = pkin(4) * t79 + t158;
t29 = -pkin(5) * t79 + t41;
t21 = pkin(4) * t75 + t141;
t15 = pkin(5) * t74 + t25;
t14 = -pkin(5) * t75 - t24;
t10 = t217 * t75 + t141;
t8 = t130 + t216;
t5 = -t117 * t17 + t120 * t13;
t1 = t120 * t2;
t18 = [qJDD(1), t228, t155, qJDD(2) - t228 - 0.2e1 * t198, 0.2e1 * t110 + t218 - t155, -t161 * pkin(1) - g(1) * (-pkin(1) * t119 + t103) - t171 + (t110 + t218) * qJ(2), t128 * t113, t128 * t114, t228 + t184 * (-t221 - t80) t84 * qJ(2) + t97 * qJD(2) - g(1) * (t116 * t119 + t103) - g(2) * (qJ(3) * t122 + t186) + t116 * t167 - qJD(3) * t229, t150, t38 * t79 + t39 * t78 - t70 * t74 - t72 * t75, t203, t144, 0, qJD(2) * t70 + t39 * t93 + t75 * t83 + t77 * t79 + t223, qJD(2) * t72 - t38 * t93 + t74 * t83 - t77 * t78 + t222, t24 * t70 + t25 * t72 - t38 * t40 - t39 * t41 - t129, -t21 * t70 - t27 * t75 - t35 * t39 - t79 * t8 - t223, -t21 * t72 - t27 * t74 + t35 * t38 + t78 * t8 - t222, t8 * t35 + t27 * t21 - t7 * t41 + t26 * t24 + t9 * t40 + t23 * t25 - g(1) * t103 - t171 + (-g(1) * t137 + g(2) * t115) * t122 + (-g(1) * (-pkin(1) + t115) - g(2) * t137) * t119, t117 * t152 + t174 * t47 (-t117 * t45 + t120 * t47) * t75 + (t201 - t117 * t12 + (-t117 * t47 - t120 * t45) * qJD(6)) * t79, -t11 * t78 + t117 * t151 + t174 * t64 + t47 * t74, -t180 * t64 * t79 + t12 * t78 + t120 * t151 - t45 * t74, t37 * t78 + t64 * t74, -g(1) * t66 - g(2) * t68 - t1 * t78 + t29 * t12 + t14 * t45 + t5 * t74 + (-t10 * t64 + t4 * t78 + t211) * t117 + (t15 * t64 - t138) * t120 + ((-t117 * t28 - t120 * t22) * t64 + t6 * t78 + t117 * t212) * qJD(6), -g(1) * t65 - g(2) * t67 + t29 * t11 + t14 * t47 - t6 * t74 + (-(qJD(6) * t28 + t10) * t64 + t211 + t169 * t78 + qJD(6) * t212) * t120 + (-(-qJD(6) * t22 + t15) * t64 + t168 * t78 + t138) * t117; 0, 0, 0, qJDD(1), -t124, -qJ(2) * t124 + t161 - t228, -t124 * t113, -t124 * t114, -t184 * qJDD(1), -qJD(1) * t97 + t167 - t228, 0, 0, 0, 0, 0, -t226, -t135, -t39 * t79 - t70 * t75 - t150, t226, t135, -qJD(1) * t27 + t129, 0, 0, 0, 0, 0, -t78 * t32 + t12 * t79 + t45 * t75 + (t117 * t163 - t120 * t74) * t64, t78 * t200 + (t117 * t74 + t120 * t163) * t64 + t152; 0, 0, 0, 0, 0, 0, t176, t175, -t184 * t124, qJD(1) * t229 + t134, 0, 0, 0, 0, 0, t233, -t149 - 0.2e1 * t195, -t219 - t220, -t233, t195 + t38, t216 + t202 - t26 * t70 + t95 + (-qJD(5) - t23) * t72 + t134, 0, 0, 0, 0, 0, t136 + t209, t117 * t230 + t207 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t219 - t220, -t149 (t72 - t173) * qJD(4) - t133, qJDD(4), -t72 * t83 + t127 + t196, qJD(4) * t204 + t70 * t83 - t232, pkin(4) * t38 - qJ(5) * t39 + (-t26 - t31) * t72 + (t23 - t231) * t70, t36 * t70 + t126 - 0.2e1 * t194 - t196, 0.2e1 * t177 - t27 * t70 + t36 * t72 + (0.2e1 * qJD(5) - t204) * qJD(4) + t232, -t9 * pkin(4) - g(3) * t156 - t7 * qJ(5) - t231 * t26 - t23 * t31 - t27 * t36 - t228 * (pkin(4) * t99 + qJ(5) * t98) -t164 * t47 + t201 (-t12 - t208) * t120 + (-t11 + t210) * t117, t142 + t207, t136 - t209, t64 * t70, qJ(5) * t12 + t125 * t117 + t224 * t120 + t187 * t45 + t5 * t70, qJ(5) * t11 - t224 * t117 + t125 * t120 + t187 * t47 - t6 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195 - t38, qJDD(4) - t206, -qJD(4) ^ 2 - t219, qJD(4) * t26 + t126 - t194, 0, 0, 0, 0, 0, -qJD(4) * t45 + t142, -qJD(4) * t47 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t11 + t210, -t12 + t208, -t37, -g(1) * t67 + g(2) * t65 - t117 * t4 - t120 * t214 - t16 * t47 + t1 + (-qJD(6) + t64) * t6, g(1) * t68 - g(2) * t66 + t16 * t45 + t5 * t64 - t169 * t120 + (-t168 + t214) * t117;];
tau_reg  = t18;
