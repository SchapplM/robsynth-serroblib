% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:17
% EndTime: 2019-03-08 23:59:27
% DurationCPUTime: 3.19s
% Computational Cost: add. (4402->319), mult. (10845->455), div. (0->0), fcn. (8166->10), ass. (0->191)
t152 = sin(qJ(3));
t153 = sin(qJ(2));
t148 = sin(pkin(6));
t220 = qJD(1) * t148;
t204 = t153 * t220;
t244 = qJD(3) * pkin(3);
t260 = -t152 * t244 + t204;
t151 = sin(qJ(4));
t155 = cos(qJ(3));
t255 = cos(qJ(4));
t167 = -t151 * t152 + t255 * t155;
t211 = qJD(3) + qJD(4);
t91 = t211 * t167;
t225 = t151 * t155;
t122 = t255 * t152 + t225;
t92 = t211 * t122;
t269 = pkin(4) * t92 - pkin(10) * t91 - t260;
t156 = cos(qJ(2));
t200 = t156 * t220;
t100 = t167 * t200;
t256 = -pkin(9) - pkin(8);
t206 = qJD(3) * t256;
t124 = t152 * t206;
t125 = t155 * t206;
t130 = t256 * t152;
t131 = t256 * t155;
t259 = t255 * t130 + t151 * t131;
t52 = t259 * qJD(4) + t255 * t124 + t151 * t125;
t268 = t52 - t100;
t197 = qJD(4) * t255;
t192 = -t256 * qJD(2) + t204;
t149 = cos(pkin(6));
t219 = qJD(1) * t149;
t97 = t152 * t219 + t192 * t155;
t86 = t151 * t97;
t96 = -t192 * t152 + t155 * t219;
t50 = t255 * t96 - t86;
t266 = -pkin(3) * t197 + t50;
t150 = sin(qJ(5));
t154 = cos(qJ(5));
t267 = t100 * t150 + t154 * t269;
t213 = qJD(5) * t154;
t143 = -pkin(3) * t155 - pkin(2);
t83 = -pkin(4) * t167 - pkin(10) * t122 + t143;
t265 = t150 * t269 + t268 * t154 + t83 * t213;
t175 = -qJ(6) * t91 - qJD(6) * t122;
t106 = t151 * t130 - t255 * t131;
t95 = t154 * t106;
t264 = t92 * pkin(5) - t150 * t52 + t175 * t154 + (-t95 + (qJ(6) * t122 - t83) * t150) * qJD(5) + t267;
t201 = t122 * t213;
t263 = -qJ(6) * t201 + (-qJD(5) * t106 + t175) * t150 + t265;
t249 = t106 * qJD(4) - t122 * t200 + t151 * t124 - t255 * t125;
t215 = qJD(4) * t151;
t218 = qJD(2) * t148;
t196 = qJD(1) * t218;
t184 = t156 * t196;
t64 = t96 * qJD(3) + t155 * t184;
t65 = -t97 * qJD(3) - t152 * t184;
t88 = t96 + t244;
t16 = t151 * t64 + t97 * t197 + t88 * t215 - t255 * t65;
t214 = qJD(5) * t150;
t47 = t255 * t88 - t86;
t40 = -t211 * pkin(4) - t47;
t262 = t16 * t154 - t40 * t214;
t198 = qJD(2) * t255;
t217 = qJD(2) * t152;
t116 = t151 * t217 - t155 * t198;
t233 = t116 * t150;
t261 = -qJ(6) * t233 + t154 * qJD(6);
t118 = -qJD(2) * t225 - t152 * t198;
t82 = -pkin(4) * t118 + pkin(10) * t116;
t70 = pkin(3) * t217 + t82;
t258 = t150 * t70 + t154 * t266;
t160 = t91 * qJD(2);
t165 = t154 * t118 - t150 * t211;
t45 = -t165 * qJD(5) + t150 * t160;
t257 = t165 ^ 2;
t254 = t154 * pkin(5);
t253 = -qJ(6) - pkin(10);
t112 = qJD(5) + t116;
t87 = t255 * t97;
t48 = t151 * t88 + t87;
t41 = t211 * pkin(10) + t48;
t108 = t143 * qJD(2) - t200;
t60 = t116 * pkin(4) + t118 * pkin(10) + t108;
t25 = -t150 * t41 + t154 * t60;
t17 = qJ(6) * t165 + t25;
t10 = pkin(5) * t112 + t17;
t252 = t10 - t17;
t251 = t150 * t82 + t154 * t47;
t141 = pkin(3) * t151 + pkin(10);
t222 = -qJ(6) - t141;
t191 = qJD(5) * t222;
t248 = t150 * t191 - t258 + t261;
t145 = t154 * qJ(6);
t181 = -t118 * pkin(5) + t116 * t145;
t67 = t154 * t70;
t247 = t154 * t191 - t181 - t67 + (-qJD(6) + t266) * t150;
t246 = t150 * t83 + t95;
t245 = qJD(2) * pkin(2);
t80 = t92 * qJD(2);
t243 = t150 * t80;
t242 = t154 * t80;
t241 = t154 * t91;
t189 = t154 * t211;
t44 = -qJD(5) * t189 - t118 * t214 - t154 * t160;
t239 = t44 * t150;
t193 = qJD(5) * t253;
t238 = t150 * t193 - t251 + t261;
t194 = -t150 * t47 + t154 * t82;
t237 = -t150 * qJD(6) + t154 * t193 - t181 - t194;
t101 = -t150 * t118 - t189;
t236 = t101 * t112;
t235 = t165 * t112;
t234 = t112 * t118;
t232 = t116 * t154;
t231 = t118 * t116;
t230 = t122 * t150;
t229 = t122 * t154;
t228 = t148 * t153;
t227 = t148 * t156;
t158 = qJD(2) ^ 2;
t226 = t148 * t158;
t157 = qJD(3) ^ 2;
t224 = t157 * t152;
t223 = t157 * t155;
t212 = qJD(2) * qJD(3);
t195 = t152 * t212;
t113 = pkin(3) * t195 + t153 * t196;
t221 = t152 ^ 2 - t155 ^ 2;
t216 = qJD(2) * t153;
t207 = t153 * t226;
t203 = t148 * t216;
t202 = t156 * t218;
t190 = t112 * t154;
t187 = t152 * t202;
t186 = t155 * t202;
t142 = -t255 * pkin(3) - pkin(4);
t26 = t150 * t60 + t154 * t41;
t185 = -t26 * t118 + t16 * t150 + t40 * t213;
t49 = t151 * t96 + t87;
t183 = pkin(3) * t215 - t49;
t182 = (t214 + t233) * pkin(5);
t18 = -qJ(6) * t101 + t26;
t179 = -t10 * t154 - t150 * t18;
t178 = t116 * t40 - t141 * t80;
t7 = pkin(5) * t45 + t16;
t174 = t25 * t118 - t262;
t114 = t149 * t155 - t152 * t228;
t115 = t149 * t152 + t155 * t228;
t72 = t151 * t114 + t255 * t115;
t58 = -t150 * t72 - t154 * t227;
t173 = t150 * t227 - t154 * t72;
t172 = t108 * t118 - t16;
t171 = t150 * t91 + t201;
t170 = -t122 * t214 + t241;
t15 = t151 * t65 + t88 * t197 - t97 * t215 + t255 * t64;
t34 = t80 * pkin(4) - pkin(10) * t160 + t113;
t169 = t154 * t15 + t150 * t34 + t60 * t213 - t41 * t214;
t168 = t255 * t114 - t151 * t115;
t166 = qJD(2) * t245;
t164 = -0.2e1 * qJD(3) * t245;
t33 = t154 * t34;
t163 = -t26 * qJD(5) - t150 * t15 + t33;
t1 = t80 * pkin(5) + t44 * qJ(6) + qJD(6) * t165 + t163;
t3 = -qJ(6) * t45 - qJD(6) * t101 + t169;
t162 = t179 * qJD(5) - t1 * t150 - t10 * t232 + t3 * t154 - t18 * t233;
t161 = t108 * t116 - t15;
t129 = pkin(10) * t154 + t145;
t128 = t253 * t150;
t120 = t141 * t154 + t145;
t119 = t222 * t150;
t98 = t101 ^ 2;
t94 = -t115 * qJD(3) - t187;
t93 = t114 * qJD(3) + t186;
t78 = t154 * t83;
t61 = -t116 ^ 2 + t118 ^ 2;
t56 = (-qJD(2) * t122 - t118) * t211;
t55 = t116 * t211 + t160;
t35 = -qJ(6) * t230 + t246;
t31 = -pkin(5) * t167 - t106 * t150 - t122 * t145 + t78;
t30 = t101 * pkin(5) + qJD(6) + t40;
t29 = t72 * qJD(4) + t151 * t93 - t255 * t94;
t28 = t168 * qJD(4) + t151 * t94 + t255 * t93;
t23 = t112 * t190 - t118 * t165 + t243;
t22 = -t112 ^ 2 * t150 - t101 * t118 + t242;
t20 = -t165 * t190 - t239;
t12 = t173 * qJD(5) - t150 * t28 + t154 * t203;
t11 = t58 * qJD(5) + t150 * t203 + t154 * t28;
t5 = (-t44 - t236) * t154 + (-t45 + t235) * t150;
t2 = [0, 0, -t207, -t156 * t226, 0, 0, 0, 0, 0, -t155 * t207 + (t94 - t187) * qJD(3), t152 * t207 + (-t93 - t186) * qJD(3), 0, 0, 0, 0, 0, -t29 * t211 + (t116 * t216 - t156 * t80) * t148, -t28 * t211 + (-t153 * t118 - t156 * t91) * t218, 0, 0, 0, 0, 0, t101 * t29 + t112 * t12 - t168 * t45 + t58 * t80, -t11 * t112 - t165 * t29 + t168 * t44 + t173 * t80, -t101 * t11 + t12 * t165 + t173 * t45 + t44 * t58, t1 * t58 + t10 * t12 + t11 * t18 - t168 * t7 - t173 * t3 + t29 * t30; 0, 0, 0, 0, 0.2e1 * t155 * t195, -0.2e1 * t221 * t212, t223, -t224, 0, -pkin(8) * t223 + t152 * t164, pkin(8) * t224 + t155 * t164, -t118 * t91 + t122 * t160, -t91 * t116 + t118 * t92 - t122 * t80 + t160 * t167, t91 * t211, -t92 * t211, 0, t108 * t92 - t113 * t167 - t260 * t116 + t143 * t80 - t249 * t211, t108 * t91 + t113 * t122 + t260 * t118 + t143 * t160 - t211 * t268, -t165 * t170 - t229 * t44 (-t101 * t154 + t150 * t165) * t91 + (t239 - t154 * t45 + (t101 * t150 + t154 * t165) * qJD(5)) * t122, t112 * t170 - t165 * t92 + t167 * t44 + t229 * t80, -t101 * t92 - t112 * t171 + t167 * t45 - t230 * t80, t112 * t92 - t167 * t80, t78 * t80 - (-t213 * t41 + t33) * t167 + t25 * t92 - t259 * t45 + t40 * t201 + (-t106 * t213 + t267) * t112 + t249 * t101 + ((-qJD(5) * t83 - t52) * t112 - t106 * t80 - (-qJD(5) * t60 - t15) * t167 + t16 * t122 + t40 * t91) * t150, -t246 * t80 + t169 * t167 - t26 * t92 + t259 * t44 + t40 * t241 + t262 * t122 + (t106 * t214 - t265) * t112 - t249 * t165, t31 * t44 - t35 * t45 + t179 * t91 + t264 * t165 - t263 * t101 + (-t1 * t154 - t150 * t3 + (t10 * t150 - t154 * t18) * qJD(5)) * t122, t3 * t35 + t1 * t31 + t7 * (pkin(5) * t230 - t259) + (pkin(5) * t171 + t249) * t30 + t263 * t18 + t264 * t10; 0, 0, 0, 0, -t152 * t158 * t155, t221 * t158, 0, 0, 0, t152 * t166, t155 * t166, -t231, t61, t55, t56, 0, t49 * t211 + (-t116 * t217 - t211 * t215) * pkin(3) + t172, t50 * t211 + (t118 * t217 - t197 * t211) * pkin(3) + t161, t20, t5, t23, t22, t234, t142 * t45 + t178 * t150 + t183 * t101 + (-t141 * t213 + t266 * t150 - t67) * t112 + t174, -t142 * t44 + t178 * t154 - t183 * t165 + (t141 * t214 + t258) * t112 + t185, -t248 * t101 + t119 * t44 - t120 * t45 + t165 * t247 + t162, t3 * t120 + t1 * t119 + t7 * (t142 - t254) + (-t87 + (pkin(3) * qJD(4) - t96) * t151 + t182) * t30 + t248 * t18 + t247 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, t61, t55, t56, 0, t211 * t48 + t172, t211 * t47 + t161, t20, t5, t23, t22, t234, -pkin(4) * t45 - t194 * t112 - t48 * t101 + t40 * t233 + (-t112 * t213 - t243) * pkin(10) + t174, pkin(4) * t44 + t251 * t112 + t48 * t165 + t40 * t232 + (t112 * t214 - t242) * pkin(10) + t185, -t238 * t101 + t128 * t44 - t129 * t45 + t165 * t237 + t162, t3 * t129 + t1 * t128 + t7 * (-pkin(4) - t254) + (t182 - t48) * t30 + t238 * t18 + t237 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t101, -t98 + t257, -t44 + t236, -t45 - t235, t80, t26 * t112 + t165 * t40 + t163, t101 * t40 + t112 * t25 - t169, pkin(5) * t44 - t252 * t101, t252 * t18 + (t165 * t30 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98 - t257, -t10 * t165 + t101 * t18 + t7;];
tauc_reg  = t2;
