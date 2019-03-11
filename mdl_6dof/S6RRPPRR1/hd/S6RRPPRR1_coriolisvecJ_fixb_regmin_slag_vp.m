% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:03
% EndTime: 2019-03-09 08:48:12
% DurationCPUTime: 2.88s
% Computational Cost: add. (3986->315), mult. (10139->414), div. (0->0), fcn. (7693->8), ass. (0->179)
t157 = qJD(2) - qJD(5);
t168 = cos(qJ(6));
t162 = sin(pkin(10));
t163 = cos(pkin(10));
t167 = sin(qJ(2));
t170 = cos(qJ(2));
t136 = t162 * t170 + t163 * t167;
t122 = t136 * qJD(2);
t105 = qJD(1) * t122;
t206 = qJD(1) * qJD(2);
t201 = t167 * t206;
t145 = t162 * t201;
t200 = t170 * t206;
t106 = t163 * t200 - t145;
t212 = qJD(1) * t170;
t202 = t163 * t212;
t213 = qJD(1) * t167;
t120 = t162 * t213 - t202;
t123 = t136 * qJD(1);
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t209 = qJD(5) * t169;
t210 = qJD(5) * t166;
t194 = -t166 * t105 - t169 * t106 - t120 * t209 + t123 * t210;
t207 = qJD(6) * t168;
t165 = sin(qJ(6));
t208 = qJD(6) * t165;
t250 = t166 * t120 + t169 * t123;
t11 = -t157 * t207 - t168 * t194 - t208 * t250;
t183 = t165 * t157 - t168 * t250;
t12 = -qJD(6) * t183 - t165 * t194;
t252 = -t169 * t120 + t166 * t123;
t258 = qJD(6) + t252;
t240 = t183 * t258;
t52 = t168 * t157 + t165 * t250;
t241 = t52 * t258;
t275 = (t11 - t241) * t168 + (-t12 + t240) * t165;
t239 = t183 * t250;
t26 = qJD(5) * t250 - t169 * t105 + t166 * t106;
t230 = t165 * t26;
t268 = t168 * t258;
t265 = t258 * t268 + t230;
t272 = t239 + t265;
t232 = t11 * t165;
t271 = t183 * t268 - t232;
t227 = t250 * t157;
t269 = t26 + t227;
t235 = -qJ(3) - pkin(7);
t143 = t235 * t167;
t140 = qJD(1) * t143;
t144 = t235 * t170;
t141 = qJD(1) * t144;
t225 = t162 * t141;
t80 = t163 * t140 + t225;
t216 = qJD(4) - t80;
t196 = -t168 * t26 + t208 * t258;
t229 = t165 * t258;
t186 = t229 * t252 + t196;
t238 = t52 * t250;
t267 = t186 - t238;
t226 = t252 * t157;
t266 = t194 + t226;
t142 = -qJD(1) * pkin(1) - pkin(2) * t212 + qJD(3);
t59 = t120 * pkin(3) - t123 * qJ(4) + t142;
t40 = -t120 * pkin(4) - t59;
t10 = pkin(5) * t252 - pkin(9) * t250 + t40;
t132 = qJD(2) * pkin(2) + t140;
t75 = t163 * t132 + t225;
t184 = qJD(4) - t75;
t245 = t123 * pkin(8);
t42 = -t245 + (-pkin(3) - pkin(4)) * qJD(2) + t184;
t246 = t120 * pkin(8);
t224 = t163 * t141;
t76 = t162 * t132 - t224;
t73 = qJD(2) * qJ(4) + t76;
t48 = t73 + t246;
t17 = t166 * t42 + t169 * t48;
t15 = -t157 * pkin(9) + t17;
t1 = t168 * t10 - t165 * t15;
t264 = t1 * t250;
t2 = t165 * t10 + t168 * t15;
t263 = t2 * t250;
t262 = t250 * pkin(5);
t261 = t258 * t250;
t260 = t250 * t252;
t259 = -t245 + t216;
t257 = t250 ^ 2 - t252 ^ 2;
t158 = qJD(2) * qJD(4);
t195 = qJD(2) * t235;
t115 = t170 * qJD(3) + t167 * t195;
t93 = t115 * qJD(1);
t116 = -t167 * qJD(3) + t170 * t195;
t94 = t116 * qJD(1);
t54 = t162 * t94 + t163 * t93;
t49 = t158 + t54;
t34 = t105 * pkin(8) + t49;
t51 = t162 * t93 - t163 * t94;
t36 = -t106 * pkin(8) + t51;
t4 = t166 * t34 - t169 * t36 + t48 * t209 + t42 * t210;
t256 = -t250 * t40 - t4;
t3 = t166 * t36 + t169 * t34 + t42 * t209 - t210 * t48;
t255 = t252 * t40 - t3;
t254 = -0.2e1 * t206;
t253 = -qJD(6) + t258;
t118 = t123 ^ 2;
t251 = -t120 ^ 2 - t118;
t61 = pkin(2) * t213 + t123 * pkin(3) + t120 * qJ(4);
t43 = -t123 * pkin(4) - t61;
t154 = -t163 * pkin(2) - pkin(3);
t149 = -pkin(4) + t154;
t152 = t162 * pkin(2) + qJ(4);
t179 = t166 * t149 + t169 * t152;
t92 = -pkin(9) + t179;
t249 = (-pkin(9) * t252 + qJD(6) * t92 - t262 + t43) * t258 - t4;
t248 = (t258 * pkin(9) + t262) * t258 + t4;
t16 = -t166 * t48 + t169 * t42;
t14 = t157 * pkin(5) - t16;
t223 = t163 * t170;
t135 = t162 * t167 - t223;
t181 = t169 * t135 - t166 * t136;
t236 = -t170 * pkin(2) - pkin(1);
t74 = t135 * pkin(3) - t136 * qJ(4) + t236;
t50 = -t135 * pkin(4) - t74;
t78 = t166 * t135 + t169 * t136;
t18 = -pkin(5) * t181 - t78 * pkin(9) + t50;
t83 = -t163 * t143 - t162 * t144;
t64 = -t136 * pkin(8) + t83;
t84 = t162 * t143 - t163 * t144;
t65 = t135 * pkin(8) + t84;
t22 = t166 * t64 + t169 * t65;
t211 = qJD(2) * t167;
t125 = qJD(2) * t223 - t162 * t211;
t31 = qJD(5) * t181 + t166 * t122 + t169 * t125;
t21 = t166 * t65 - t169 * t64;
t62 = t162 * t115 - t163 * t116;
t44 = -t125 * pkin(8) + t62;
t63 = t163 * t115 + t162 * t116;
t45 = t122 * pkin(8) + t63;
t8 = -qJD(5) * t21 + t166 * t44 + t169 * t45;
t247 = t14 * t31 - (qJD(6) * t18 + t8) * t258 + (qJD(6) * t10 + t3) * t181 - t22 * t26 + t4 * t78;
t244 = t14 * t78;
t243 = t18 * t26;
t242 = t51 * t83;
t237 = t78 * t26;
t180 = t169 * t149 - t166 * t152;
t79 = t162 * t140 - t224;
t56 = t79 + t246;
t234 = -qJD(5) * t180 + t166 * t56 - t259 * t169;
t233 = qJD(5) * t179 + t259 * t166 + t169 * t56;
t172 = qJD(1) ^ 2;
t219 = t170 * t172;
t171 = qJD(2) ^ 2;
t218 = t171 * t167;
t217 = t171 * t170;
t214 = t167 ^ 2 - t170 ^ 2;
t205 = pkin(2) * t211;
t204 = t78 * t208;
t203 = t258 * t207;
t189 = pkin(1) * t254;
t188 = t157 * t258;
t187 = t157 ^ 2;
t150 = pkin(2) * t201;
t37 = t105 * pkin(3) - t106 * qJ(4) - t123 * qJD(4) + t150;
t185 = t258 * t31 + t237;
t178 = t59 * t123 + t51;
t24 = -t105 * pkin(4) - t37;
t47 = t122 * pkin(3) - t125 * qJ(4) - t136 * qJD(4) + t205;
t176 = -pkin(9) * t26 + (t14 + t16) * t258;
t33 = -t122 * pkin(4) - t47;
t175 = -t84 * t105 + t83 * t106 - t63 * t120 + t62 * t123 + t51 * t136;
t174 = -t92 * t26 + (-t14 + t234) * t258;
t91 = pkin(5) - t180;
t67 = -qJD(2) * pkin(3) + t184;
t32 = qJD(5) * t78 - t169 * t122 + t166 * t125;
t9 = qJD(5) * t22 + t166 * t45 - t169 * t44;
t7 = t32 * pkin(5) - t31 * pkin(9) + t33;
t6 = t26 * pkin(5) + pkin(9) * t194 + t24;
t5 = t168 * t6;
t13 = [0, 0, 0, 0.2e1 * t167 * t200, t214 * t254, t217, -t218, 0, -pkin(7) * t217 + t167 * t189, pkin(7) * t218 + t170 * t189, -t76 * t122 - t75 * t125 - t54 * t135 + t175, t242 + t54 * t84 - t75 * t62 + t76 * t63 + (qJD(1) * t236 + t142) * t205, -t62 * qJD(2) + t74 * t105 + t47 * t120 + t59 * t122 + t37 * t135, -t73 * t122 + t67 * t125 - t49 * t135 + t175, t63 * qJD(2) - t74 * t106 - t47 * t123 - t59 * t125 - t37 * t136, t37 * t74 + t59 * t47 + t49 * t84 + t67 * t62 + t73 * t63 + t242, -t194 * t78 + t250 * t31, -t181 * t194 - t250 * t32 - t252 * t31 - t237, -t31 * t157, t32 * t157, 0, t9 * t157 - t181 * t24 + t252 * t33 + t50 * t26 + t40 * t32, t8 * t157 - t194 * t50 + t24 * t78 + t250 * t33 + t40 * t31, t183 * t204 + (t11 * t78 - t183 * t31) * t168 (t165 * t183 - t168 * t52) * t31 + (-t232 - t12 * t168 + (t165 * t52 + t168 * t183) * qJD(6)) * t78, -t11 * t181 + t168 * t185 - t183 * t32 - t204 * t258, t12 * t181 - t165 * t185 - t203 * t78 - t52 * t32, -t181 * t26 + t258 * t32, t1 * t32 + t21 * t12 - t5 * t181 + t9 * t52 + (t243 + t7 * t258 + (t15 * t181 - t22 * t258 + t244) * qJD(6)) * t168 + t247 * t165, t21 * t11 - t2 * t32 - t9 * t183 + (-(-qJD(6) * t22 + t7) * t258 - t243 + (-qJD(6) * t15 + t6) * t181 - qJD(6) * t244) * t165 + t247 * t168; 0, 0, 0, -t167 * t219, t214 * t172, 0, 0, 0, t172 * pkin(1) * t167, pkin(1) * t219 (t76 - t79) * t123 + (-t75 + t80) * t120 + (-t105 * t162 - t106 * t163) * pkin(2), t75 * t79 - t76 * t80 + (-t142 * t213 + t162 * t54 - t163 * t51) * pkin(2), t79 * qJD(2) - t61 * t120 - t178, -t152 * t105 + t154 * t106 + (t73 - t79) * t123 + (t67 - t216) * t120, -t80 * qJD(2) - t59 * t120 + t61 * t123 + 0.2e1 * t158 + t54, t49 * t152 + t51 * t154 + t216 * t73 - t59 * t61 - t67 * t79, -t260, -t257, t266, t269, 0, t233 * t157 - t252 * t43 - t256, -t234 * t157 - t250 * t43 - t255, t271, -t275, -t272, t267, t261, t91 * t12 + t174 * t165 - t249 * t168 + t233 * t52 + t264, t91 * t11 + t249 * t165 + t174 * t168 - t183 * t233 - t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t76 * t120 + t75 * t123 + t150, 0.2e1 * t123 * qJD(2), t251, t145 + (t120 - t202) * qJD(2), t73 * t120 - t67 * t123 + t37, 0, 0, 0, 0, 0, -t26 + t227, t194 - t226, 0, 0, 0, 0, 0, t186 + t238, -t239 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t120, -t145 + (t120 + t202) * qJD(2), -t118 - t171, -t73 * qJD(2) + t178, 0, 0, 0, 0, 0, -t123 * t252 - t166 * t187, -t123 * t250 - t169 * t187, 0, 0, 0, 0, 0, -t123 * t268 + (t165 * t188 - t12) * t169 + (-t157 * t52 - t203 - t230) * t166, t123 * t229 + (t168 * t188 - t11) * t169 + (t157 * t183 + t196) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t257, -t266, -t269, 0, -t17 * t157 + t256, -t16 * t157 + t255, -t271, t275, t272, -t267, -t261, -pkin(5) * t12 + t176 * t165 - t248 * t168 - t17 * t52 - t264, -pkin(5) * t11 + t248 * t165 + t176 * t168 + t17 * t183 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 * t52, t183 ^ 2 - t52 ^ 2, t11 + t241, -t12 - t240, t26, t14 * t183 - t165 * t3 + t253 * t2 + t5, t253 * t1 + t14 * t52 - t165 * t6 - t168 * t3;];
tauc_reg  = t13;
