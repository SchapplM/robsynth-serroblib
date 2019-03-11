% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:01
% EndTime: 2019-03-08 23:14:09
% DurationCPUTime: 2.59s
% Computational Cost: add. (3137->305), mult. (7877->418), div. (0->0), fcn. (5943->10), ass. (0->185)
t139 = sin(qJ(4));
t140 = sin(qJ(3));
t143 = cos(qJ(3));
t245 = cos(qJ(4));
t107 = t139 * t143 + t245 * t140;
t100 = t107 * qJD(2);
t258 = qJD(6) + t100;
t132 = qJD(3) + qJD(4);
t138 = sin(qJ(6));
t142 = cos(qJ(6));
t199 = t245 * t143;
t182 = qJD(2) * t199;
t210 = qJD(2) * t140;
t194 = t139 * t210;
t98 = -t182 + t194;
t84 = t138 * t132 - t142 * t98;
t260 = t258 * t84;
t185 = t138 * t258;
t259 = qJD(6) - t258;
t247 = -pkin(9) - pkin(8);
t141 = sin(qJ(2));
t136 = sin(pkin(6));
t213 = qJD(1) * t136;
t196 = t141 * t213;
t236 = qJD(3) * pkin(3);
t155 = t140 * t236 - t196;
t200 = qJD(3) * t247;
t109 = t140 * t200;
t110 = t143 * t200;
t113 = t247 * t140;
t114 = t247 * t143;
t192 = t245 * qJD(4);
t144 = cos(qJ(2));
t195 = t144 * t213;
t208 = qJD(4) * t139;
t220 = t139 * t140;
t253 = t199 - t220;
t239 = -t245 * t109 - t139 * t110 - t113 * t192 - t114 * t208 + t253 * t195;
t159 = t139 * t113 - t245 * t114;
t238 = t159 * qJD(4) - t107 * t195 + t139 * t109 - t245 * t110;
t188 = -t247 * qJD(2) + t196;
t137 = cos(pkin(6));
t212 = qJD(1) * t137;
t81 = t140 * t212 + t188 * t143;
t75 = t139 * t81;
t80 = -t188 * t140 + t143 * t212;
t38 = t245 * t80 - t75;
t224 = pkin(3) * t192 + qJD(5) - t38;
t76 = t245 * t81;
t37 = t139 * t80 + t76;
t179 = pkin(3) * t208 - t37;
t257 = t238 * t132;
t256 = t239 * t132;
t77 = t80 + t236;
t35 = -t245 * t77 + t75;
t216 = qJD(5) + t35;
t174 = t132 * t220;
t78 = -qJD(3) * t199 - t143 * t192 + t174;
t254 = -t78 * qJ(5) + t107 * qJD(5) - t155;
t223 = t136 * t141;
t161 = -t137 * t143 + t140 * t223;
t252 = qJD(3) * t161;
t211 = qJD(2) * t136;
t197 = t144 * t211;
t183 = t143 * t197;
t148 = t183 - t252;
t97 = t137 * t140 + t143 * t223;
t152 = t97 * qJD(3);
t184 = t140 * t197;
t149 = t152 + t184;
t18 = t139 * t149 - t245 * t148 + t161 * t192 + t97 * t208;
t209 = qJD(2) * t141;
t215 = t132 * t182;
t65 = qJD(2) * t174 - t215;
t251 = t136 * (-t100 * t209 - t144 * t65) - t18 * t132;
t59 = -t139 * t161 + t245 * t97;
t19 = t59 * qJD(4) + t139 * t148 + t245 * t149;
t79 = t132 * t107;
t66 = t79 * qJD(2);
t250 = t136 * (-t144 * t66 + t98 * t209) - t19 * t132;
t249 = t100 ^ 2;
t248 = pkin(4) + pkin(10);
t246 = t98 * pkin(5);
t244 = t100 * pkin(5);
t128 = -t143 * pkin(3) - pkin(2);
t171 = -t107 * qJ(5) + t128;
t47 = -t248 * t253 + t171;
t243 = t47 * t65;
t242 = t258 * t98;
t241 = -t79 * pkin(5) - t239;
t240 = -t78 * pkin(5) + t238;
t36 = t139 * t77 + t76;
t237 = qJD(2) * pkin(2);
t235 = t100 * t98;
t234 = t138 * t65;
t86 = t142 * t132 + t138 * t98;
t233 = t138 * t86;
t61 = t142 * t65;
t232 = t142 * t258;
t31 = -t132 * qJ(5) - t36;
t24 = -t31 - t246;
t229 = t24 * t253;
t206 = qJD(6) * t142;
t207 = qJD(6) * t138;
t33 = -t132 * t207 + t138 * t66 + t98 * t206;
t228 = t33 * t142;
t92 = t128 * qJD(2) - t195;
t150 = -t100 * qJ(5) + t92;
t45 = t98 * pkin(4) + t150;
t227 = t45 * t100;
t226 = t92 * t100;
t225 = t244 + t224;
t222 = t136 * t144;
t147 = qJD(2) ^ 2;
t221 = t136 * t147;
t146 = qJD(3) ^ 2;
t219 = t146 * t140;
t218 = t146 * t143;
t217 = t244 + t216;
t204 = qJD(2) * qJD(3);
t190 = t140 * t204;
t191 = qJD(1) * t211;
t95 = pkin(3) * t190 + t141 * t191;
t214 = t140 ^ 2 - t143 ^ 2;
t181 = t144 * t191;
t53 = t80 * qJD(3) + t143 * t181;
t54 = -t81 * qJD(3) - t140 * t181;
t189 = -t139 * t54 - t77 * t192 + t81 * t208 - t245 * t53;
t7 = -t132 * qJD(5) + t189;
t4 = -t66 * pkin(5) - t7;
t20 = -t248 * t132 + t217;
t32 = t248 * t98 + t150;
t9 = t138 * t20 + t142 * t32;
t202 = t4 * t142 - t9 * t98;
t201 = t141 * t221;
t198 = t136 * t209;
t10 = t139 * t53 + t81 * t192 + t77 * t208 - t245 * t54;
t187 = t258 * t24;
t127 = -t245 * pkin(3) - pkin(4);
t124 = -pkin(10) + t127;
t67 = t100 * pkin(4) + t98 * qJ(5);
t55 = pkin(3) * t210 + t67;
t93 = t100 * pkin(10);
t186 = -qJD(6) * t124 + t55 + t93;
t180 = t246 + t179;
t178 = -t248 * t79 + t254;
t177 = -t79 * pkin(4) + t254;
t173 = t138 * t32 - t142 * t20;
t176 = t4 * t138 - t173 * t98 + (t100 * t142 + t206) * t24;
t87 = -t245 * t113 - t139 * t114;
t170 = t92 * t98 + t189;
t169 = t10 + t227;
t168 = t36 * t132 - t10;
t167 = -t185 * t258 - t61;
t165 = -t207 * t258 - t61;
t58 = t139 * t97 + t161 * t245;
t164 = -t138 * t58 + t142 * t222;
t163 = t138 * t222 + t142 * t58;
t158 = t237 * qJD(2);
t157 = -t45 * t98 - t7;
t56 = t107 * pkin(5) + t87;
t156 = -t24 * t79 + t253 * t4 - t56 * t65;
t154 = t65 * qJ(5) - t100 * qJD(5) + t95;
t153 = -t232 * t258 + t234;
t151 = -0.2e1 * qJD(3) * t237;
t42 = t215 + (-t194 + t98) * t132;
t125 = t139 * pkin(3) + qJ(5);
t70 = -pkin(4) * t253 + t171;
t62 = t142 * t66;
t57 = pkin(5) * t253 + t159;
t48 = -t98 ^ 2 + t249;
t46 = t65 * t107;
t44 = t67 + t93;
t34 = t86 * qJD(6) - t62;
t30 = -t132 * pkin(4) + t216;
t27 = t36 - t246;
t16 = t66 * pkin(4) + t154;
t14 = -t84 * t98 + t153;
t13 = t86 * t98 + t167;
t12 = t248 * t66 + t154;
t11 = -t86 * t185 + t228;
t6 = -pkin(5) * t65 + t10;
t5 = t142 * t6;
t1 = (-t258 * t86 - t34) * t142 + (-t33 + t260) * t138;
t2 = [0, 0, -t201, -t144 * t221, 0, 0, 0, 0, 0, -t143 * t201 + (-t152 - 0.2e1 * t184) * qJD(3), t140 * t201 + (-0.2e1 * t183 + t252) * qJD(3), 0, 0, 0, 0, 0, t250, -t251, t19 * t100 + t18 * t98 - t58 * t65 - t59 * t66, -t250, t251, t10 * t58 + t31 * t18 + t30 * t19 - t7 * t59 + (-t144 * t16 + t209 * t45) * t136, 0, 0, 0, 0, 0 (qJD(6) * t164 - t138 * t198 + t142 * t19) * t258 - t163 * t65 - t18 * t84 + t59 * t34 -(qJD(6) * t163 + t138 * t19 + t142 * t198) * t258 - t164 * t65 - t18 * t86 + t59 * t33; 0, 0, 0, 0, 0.2e1 * t143 * t190, -0.2e1 * t214 * t204, t218, -t219, 0, -pkin(8) * t218 + t140 * t151, pkin(8) * t219 + t143 * t151, -t100 * t78 - t46, -t100 * t79 - t107 * t66 - t253 * t65 + t78 * t98, -t78 * t132, -t79 * t132, 0, t128 * t66 + t155 * t98 - t253 * t95 + t92 * t79 - t257, t155 * t100 + t95 * t107 - t128 * t65 - t92 * t78 + t256, t10 * t107 + t238 * t100 - t159 * t66 + t239 * t98 - t253 * t7 - t30 * t78 + t31 * t79 - t87 * t65, t16 * t253 + t177 * t98 - t45 * t79 - t70 * t66 + t257, t177 * t100 - t16 * t107 + t45 * t78 + t70 * t65 - t256, t10 * t87 - t159 * t7 + t16 * t70 - t177 * t45 + t238 * t30 + t239 * t31, t79 * t233 - (t33 * t138 + t206 * t86) * t253 (-t138 * t84 + t142 * t86) * t79 - (-t138 * t34 + t228 + (-t142 * t84 - t233) * qJD(6)) * t253, t79 * t185 + t33 * t107 - t86 * t78 - (t206 * t258 - t234) * t253, -t34 * t107 - t165 * t253 + t79 * t232 + t84 * t78, -t258 * t78 - t46, t5 * t107 + t57 * t34 + t173 * t78 + t241 * t84 + (-t12 * t107 + t178 * t258 + t243) * t138 + (t240 * t258 + t156) * t142 + ((-t138 * t56 - t142 * t47) * t258 - t9 * t107 - t138 * t229) * qJD(6), t57 * t33 + t9 * t78 + t241 * t86 + (t243 - (qJD(6) * t20 + t12) * t107 - qJD(6) * t229 + (-qJD(6) * t56 + t178) * t258) * t142 + (-(-qJD(6) * t32 + t6) * t107 + (qJD(6) * t47 - t240) * t258 - t156) * t138; 0, 0, 0, 0, -t140 * t147 * t143, t214 * t147, 0, 0, 0, t140 * t158, t143 * t158, t235, t48, t42, 0, 0, -t226 + t37 * t132 + (-t132 * t208 - t98 * t210) * pkin(3) - t10, t38 * t132 + (-t100 * t210 - t132 * t192) * pkin(3) + t170, -t125 * t66 - t127 * t65 + (t30 - t224) * t98 + (t179 - t31) * t100, t132 * t179 + t55 * t98 + t169, t55 * t100 + t224 * t132 + t157, t10 * t127 - t7 * t125 + t179 * t30 - t224 * t31 - t45 * t55, t11, t1, t13, t14, t242, -t124 * t61 + t125 * t34 + t225 * t84 + (t138 * t186 + t142 * t180) * t258 + t176, t125 * t33 + t225 * t86 + t186 * t232 + (t124 * t65 - t180 * t258 - t187) * t138 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t48, t42, 0, 0, t168 - t226, -t35 * t132 + t170, pkin(4) * t65 - qJ(5) * t66 + (-t31 - t36) * t100 + (t30 - t216) * t98, t67 * t98 - t168 + t227, t67 * t100 + t132 * t216 + t157, -pkin(4) * t10 - qJ(5) * t7 - t216 * t31 - t30 * t36 - t45 * t67, t11, t1, t13, t14, t242, qJ(5) * t34 - (-t138 * t44 + t142 * t27) * t258 + t217 * t84 - t165 * t248 + t176, qJ(5) * t33 + t217 * t86 + (qJD(6) * t248 + t44) * t232 + (-t248 * t65 + t258 * t27 - t187) * t138 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t235, -t132 ^ 2 - t249, t31 * t132 + t169, 0, 0, 0, 0, 0, -t132 * t84 + t167, -t132 * t86 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t84, -t84 ^ 2 + t86 ^ 2, t33 + t260, -t259 * t86 + t62, -t65, -t138 * t12 - t24 * t86 - t259 * t9 + t5, -t142 * t12 - t138 * t6 + t259 * t173 + t24 * t84;];
tauc_reg  = t2;
