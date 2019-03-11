% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:39
% EndTime: 2019-03-08 20:43:46
% DurationCPUTime: 2.86s
% Computational Cost: add. (2918->358), mult. (6014->499), div. (0->0), fcn. (4733->14), ass. (0->200)
t143 = sin(qJ(5));
t144 = sin(qJ(4));
t147 = cos(qJ(5));
t148 = cos(qJ(4));
t94 = t143 * t148 + t144 * t147;
t88 = t94 * qJD(2);
t268 = qJD(6) + t88;
t134 = qJD(4) + qJD(5);
t142 = sin(qJ(6));
t146 = cos(qJ(6));
t221 = qJD(2) * t148;
t222 = qJD(2) * t144;
t87 = t143 * t222 - t147 * t221;
t64 = -t146 * t134 - t142 * t87;
t269 = t268 * t64;
t262 = t268 - qJD(6);
t137 = qJ(4) + qJ(5);
t130 = sin(t137);
t131 = cos(t137);
t141 = cos(pkin(6));
t139 = sin(pkin(6));
t149 = cos(qJ(2));
t233 = t139 * t149;
t140 = cos(pkin(11));
t236 = t139 * t140;
t138 = sin(pkin(11));
t237 = t138 * t139;
t145 = sin(qJ(2));
t231 = t141 * t149;
t80 = t138 * t145 - t140 * t231;
t82 = t138 * t231 + t140 * t145;
t166 = -g(3) * (-t130 * t141 - t131 * t233) - g(2) * (t130 * t236 + t131 * t80) - g(1) * (-t130 * t237 + t131 * t82);
t133 = qJDD(4) + qJDD(5);
t208 = t141 * qJDD(1);
t150 = -pkin(2) - pkin(8);
t234 = t139 * t145;
t200 = qJD(2) * t234;
t101 = qJD(1) * t200;
t212 = qJDD(1) * t139;
t193 = t149 * t212;
t167 = qJDD(3) + t101 - t193;
t67 = t150 * qJDD(2) + t167;
t176 = -t144 * t208 + t148 * t67;
t225 = qJD(1) * t139;
t197 = t149 * t225;
t175 = qJD(3) - t197;
t86 = t150 * qJD(2) + t175;
t192 = pkin(9) * qJD(2) - t86;
t224 = qJD(1) * t141;
t201 = t148 * t224;
t206 = t148 * qJDD(2);
t19 = -pkin(9) * t206 + qJDD(4) * pkin(4) + (t192 * t144 - t201) * qJD(4) + t176;
t202 = t144 * t224;
t22 = t148 * t208 + (-pkin(9) * qJDD(2) + t67) * t144 + (-t192 * t148 - t202) * qJD(4);
t54 = -pkin(9) * t222 + t144 * t86 + t201;
t239 = t147 * t54;
t53 = -pkin(9) * t221 + t148 * t86 - t202;
t45 = qJD(4) * pkin(4) + t53;
t24 = t143 * t45 + t239;
t260 = t24 * qJD(5) + t143 * t22 - t147 * t19;
t4 = -pkin(5) * t133 + t260;
t164 = t166 - t4;
t213 = qJD(2) * qJD(4);
t267 = -t144 * t213 + t206;
t232 = t141 * t145;
t81 = t138 * t149 + t140 * t232;
t83 = -t138 * t232 + t140 * t149;
t162 = g(1) * t83 + g(2) * t81 + g(3) * t234;
t160 = t162 * t130;
t217 = qJD(5) * t143;
t198 = t144 * t217;
t170 = -qJD(2) * t198 + t143 * t267;
t183 = t134 * t148;
t207 = t144 * qJDD(2);
t36 = (qJD(2) * t183 + t207) * t147 + t170;
t33 = qJDD(6) + t36;
t121 = t144 * pkin(4) + qJ(3);
t93 = t143 * t144 - t147 * t148;
t47 = pkin(5) * t94 + pkin(10) * t93 + t121;
t266 = t47 * t33 - t160;
t265 = t134 * t88;
t102 = t146 * t234;
t168 = -t141 * t148 + t144 * t233;
t84 = -t141 * t144 - t148 * t233;
t40 = t143 * t84 - t147 * t168;
t264 = -t142 * t40 + t102;
t218 = qJD(4) * t148;
t177 = -pkin(4) * t218 - t175;
t181 = g(1) * t82 + g(2) * t80;
t263 = -g(3) * t233 + t181;
t241 = t143 * t54;
t23 = t147 * t45 - t241;
t20 = -pkin(5) * t134 - t23;
t203 = t145 * t225;
t249 = pkin(9) - t150;
t97 = t249 * t144;
t98 = t249 * t148;
t60 = -t143 * t97 + t147 * t98;
t219 = qJD(4) * t144;
t90 = t249 * t219;
t91 = qJD(4) * t98;
t248 = t60 * qJD(5) - t143 * t90 + t147 * t91 + t94 * t203;
t259 = (qJD(5) * t45 + t22) * t147 + t143 * t19 - t54 * t217;
t3 = pkin(10) * t133 + t259;
t223 = qJD(2) * qJ(3);
t95 = t203 + t223;
t75 = pkin(4) * t222 + t95;
t34 = pkin(5) * t88 + pkin(10) * t87 + t75;
t216 = qJD(5) * t147;
t55 = -t143 * t218 - t144 * t216 - t147 * t219 - t148 * t217;
t61 = -t143 * t98 - t147 * t97;
t261 = -(qJD(6) * t34 + t3) * t94 + t20 * t55 + (-qJD(6) * t47 + t248) * t268 - t4 * t93 - t61 * t33 + t263;
t258 = qJD(4) * (t203 - t95 - t223) - qJDD(4) * t150;
t254 = t20 * t88;
t253 = t20 * t93;
t251 = t268 * t87;
t250 = t87 * t88;
t247 = t61 * qJD(5) - t143 * t91 - t147 * t90 - t93 * t203;
t246 = -t93 * t133 + t55 * t134;
t214 = qJD(6) * t146;
t215 = qJD(6) * t142;
t35 = -t143 * t207 + t147 * t206 - t265;
t14 = t142 * t133 + t134 * t214 + t146 * t35 + t87 * t215;
t245 = t14 * t142;
t244 = t142 * t33;
t242 = t142 * t268;
t240 = t146 * t33;
t238 = qJDD(2) * pkin(2);
t235 = t139 * t144;
t152 = qJD(2) ^ 2;
t230 = t149 * t152;
t229 = t95 * qJD(2);
t136 = t148 ^ 2;
t227 = t144 ^ 2 - t136;
t151 = qJD(4) ^ 2;
t226 = -t151 - t152;
t220 = qJD(2) * t149;
t211 = qJDD(2) * qJ(3);
t210 = qJDD(4) * t144;
t205 = t93 * t215;
t204 = t142 * t234;
t199 = t139 * t220;
t21 = pkin(10) * t134 + t24;
t173 = t142 * t21 - t146 * t34;
t196 = -t173 * t87 + t20 * t215;
t194 = t148 * t213;
t187 = -t146 * t133 + t142 * t35;
t186 = t146 * t268;
t123 = pkin(4) * t143 + pkin(10);
t46 = -pkin(5) * t87 + pkin(10) * t88;
t185 = pkin(4) * t221 + qJD(6) * t123 + t46;
t184 = qJD(6) * t94 + qJD(2);
t25 = t143 * t53 + t239;
t182 = pkin(4) * t217 - t25;
t56 = -t143 * t219 + t147 * t183 - t198;
t180 = pkin(5) * t56 - pkin(10) * t55 - t177;
t66 = t134 * t142 - t146 * t87;
t179 = -t14 * t93 + t55 * t66;
t178 = -t268 * t55 + t33 * t93;
t174 = -t133 * t94 - t134 * t56;
t11 = t142 * t34 + t146 * t21;
t39 = -t143 * t168 - t147 * t84;
t172 = (-qJD(2) * pkin(2) + t175) * t145 + t95 * t149;
t171 = qJDD(2) * t145 + t230;
t169 = t146 * t40 + t204;
t165 = -t11 * t87 - t142 * t164 + t20 * t214;
t106 = t145 * t212;
t161 = -t106 + t162;
t159 = t193 + t263;
t158 = qJDD(3) - t159;
t68 = t211 + t106 + (qJD(3) + t197) * qJD(2);
t26 = t147 * t53 - t241;
t157 = -t123 * t33 + t254 + (-pkin(4) * t216 + t26) * t268;
t44 = t68 + (t194 + t207) * pkin(4);
t50 = t130 * t82 + t131 * t237;
t52 = -t80 * t130 + t131 * t236;
t73 = -t130 * t233 + t131 * t141;
t156 = g(1) * t50 - g(2) * t52 + g(3) * t73 + t75 * t88 - t259;
t155 = t75 * t87 + t166 - t260;
t153 = t175 * qJD(2) - t150 * t151 - t162 + t211 + t68;
t129 = qJDD(4) * t148;
t124 = -pkin(4) * t147 - pkin(5);
t79 = t171 * t139;
t78 = (-qJDD(2) * t149 + t145 * t152) * t139;
t70 = t167 - t238;
t58 = qJD(4) * t168 + t148 * t200;
t57 = qJD(4) * t84 + t144 * t200;
t37 = t87 ^ 2 - t88 ^ 2;
t29 = -t134 * t87 + (-t134 * t221 - t207) * t147 - t170;
t28 = t35 + t265;
t15 = qJD(6) * t66 + t187;
t13 = t40 * qJD(5) + t143 * t57 - t147 * t58;
t12 = -t39 * qJD(5) + t143 * t58 + t147 * t57;
t9 = pkin(5) * t36 - pkin(10) * t35 + t44;
t8 = t146 * t9;
t7 = t186 * t268 + t66 * t87 + t244;
t6 = -t242 * t268 - t64 * t87 + t240;
t5 = t66 * t186 + t245;
t1 = (t14 - t269) * t146 + (-t268 * t66 - t15) * t142;
t2 = [qJDD(1) - g(3), 0, -t78, -t79, t78, t79, qJDD(1) * t141 ^ 2 - g(3) + (t172 * qJD(2) + t145 * t68 - t149 * t70) * t139, 0, 0, 0, 0, 0, qJD(4) * t58 + qJDD(4) * t84 + (t144 * t171 + t145 * t194) * t139, -qJD(4) * t57 + qJDD(4) * t168 + (t145 * t267 + t148 * t230) * t139, 0, 0, 0, 0, 0, -t13 * t134 - t133 * t39 + (t145 * t36 + t88 * t220) * t139, -t12 * t134 - t133 * t40 + (t145 * t35 - t87 * t220) * t139, 0, 0, 0, 0, 0 (-qJD(6) * t169 - t12 * t142 + t146 * t199) * t268 + t264 * t33 + t13 * t64 + t39 * t15 -(qJD(6) * t264 + t146 * t12 + t142 * t199) * t268 - t169 * t33 + t13 * t66 + t39 * t14; 0, qJDD(2), t159, t161, t158 - 0.2e1 * t238, 0.2e1 * qJD(2) * qJD(3) - t161 + 0.2e1 * t211, t68 * qJ(3) + t95 * qJD(3) - t70 * pkin(2) - g(1) * (-pkin(2) * t82 + qJ(3) * t83) - g(2) * (-pkin(2) * t80 + qJ(3) * t81) + (-g(3) * (pkin(2) * t149 + qJ(3) * t145) - t172 * qJD(1)) * t139, qJDD(2) * t136 - 0.2e1 * t144 * t194, -0.2e1 * t144 * t206 + 0.2e1 * t227 * t213, -t144 * t151 + t129, -t148 * t151 - t210, 0, t153 * t144 - t148 * t258, t144 * t258 + t153 * t148, -t35 * t93 - t55 * t87, -t35 * t94 + t36 * t93 - t55 * t88 + t56 * t87, t246, t174, 0, t121 * t36 - t133 * t60 - t247 * t134 - t177 * t88 + t44 * t94 + t56 * t75 - t160, t121 * t35 - t162 * t131 - t133 * t61 + t248 * t134 + t177 * t87 - t44 * t93 + t55 * t75, t179 * t146 + t66 * t205 (-t142 * t66 - t146 * t64) * t55 + (t245 + t146 * t15 + (-t142 * t64 + t146 * t66) * qJD(6)) * t93, t14 * t94 - t178 * t146 + t205 * t268 + t56 * t66, t214 * t268 * t93 + t178 * t142 - t15 * t94 - t56 * t64, t268 * t56 + t33 * t94, -t173 * t56 + t60 * t15 + t8 * t94 + t247 * t64 + (t180 * t268 + (-t21 * t94 - t268 * t61 - t253) * qJD(6) + t266) * t146 + t261 * t142, -t11 * t56 + t60 * t14 + t247 * t66 + (-(-qJD(6) * t21 + t9) * t94 + qJD(6) * t253 + (qJD(6) * t61 - t180) * t268 - t266) * t142 + t261 * t146; 0, 0, 0, 0, qJDD(2), -t152, t101 + t158 - t229 - t238, 0, 0, 0, 0, 0, t226 * t144 + t129, t226 * t148 - t210, 0, 0, 0, 0, 0, -qJD(2) * t88 + t246, qJD(2) * t87 + t174, 0, 0, 0, 0, 0, -t94 * t244 + t15 * t93 - t55 * t64 + (-t142 * t56 - t184 * t146) * t268, -t94 * t240 + (t142 * t184 - t146 * t56) * t268 - t179; 0, 0, 0, 0, 0, 0, 0, t148 * t152 * t144, -t227 * t152, t206, -t207, qJDD(4), -t95 * t221 - g(1) * (-t138 * t235 + t148 * t82) - g(2) * (t140 * t235 + t148 * t80) - g(3) * t84 + t176, -g(3) * t168 + (-t208 + (g(1) * t138 - g(2) * t140) * t139) * t148 + (t181 - t67 + t229) * t144, -t250, t37, t28, t29, t133, t134 * t25 + (t133 * t147 - t134 * t217 - t88 * t221) * pkin(4) + t155, t134 * t26 + (-t133 * t143 - t134 * t216 + t87 * t221) * pkin(4) + t156, t5, t1, t7, t6, t251, t124 * t15 + t182 * t64 + t157 * t142 + (-t185 * t268 + t164) * t146 + t196, t124 * t14 + t146 * t157 + t182 * t66 + t185 * t242 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t37, t28, t29, t133, t134 * t24 + t155, t134 * t23 + t156, t5, t1, t7, t6, t251, -pkin(5) * t15 - t24 * t64 + (-pkin(10) * t33 + t23 * t268 + t254) * t142 + ((-pkin(10) * qJD(6) - t46) * t268 + t164) * t146 + t196, -pkin(5) * t14 + (t142 * t46 + t146 * t23) * t268 - t24 * t66 + t146 * t254 + (t215 * t268 - t240) * pkin(10) + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, t14 + t269, t262 * t66 - t187, t33, -t142 * t3 + t8 - t20 * t66 - g(1) * (-t142 * t50 + t146 * t83) - g(2) * (t142 * t52 + t146 * t81) - g(3) * (-t142 * t73 + t102) + t262 * t11, -t146 * t3 - t142 * t9 + t20 * t64 - g(1) * (-t142 * t83 - t146 * t50) - g(2) * (-t142 * t81 + t146 * t52) - g(3) * (-t146 * t73 - t204) - t262 * t173;];
tau_reg  = t2;
