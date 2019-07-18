% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:18:42
% EndTime: 2019-07-18 17:18:57
% DurationCPUTime: 3.72s
% Computational Cost: add. (6473->383), mult. (16296->556), div. (0->0), fcn. (12518->8), ass. (0->186)
t142 = qJD(2) + qJD(3);
t147 = sin(qJ(3));
t151 = cos(qJ(2));
t257 = cos(qJ(3));
t195 = qJD(1) * t257;
t148 = sin(qJ(2));
t218 = qJD(1) * t148;
t114 = t147 * t218 - t151 * t195;
t146 = sin(qJ(4));
t149 = cos(qJ(5));
t223 = t146 * t149;
t145 = sin(qJ(5));
t150 = cos(qJ(4));
t224 = t145 * t150;
t121 = t223 + t224;
t208 = qJD(4) + qJD(5);
t265 = t208 * t121;
t267 = t121 * t114 + t265;
t210 = qJD(5) * t149;
t213 = qJD(4) * t150;
t225 = t145 * t146;
t258 = t149 * t150 - t225;
t266 = t258 * t114 + t149 * t213 + t150 * t210 - t208 * t225;
t122 = t147 * t151 + t257 * t148;
t116 = t122 * qJD(1);
t97 = t116 * t146 - t150 * t142;
t99 = t116 * t150 + t142 * t146;
t175 = t145 * t97 - t149 * t99;
t60 = t145 * t99 + t149 * t97;
t254 = t60 * t175;
t264 = t175 ^ 2 - t60 ^ 2;
t194 = t257 * qJD(2);
t185 = pkin(1) * t194;
t171 = qJD(3) * t185;
t222 = t147 * t148;
t170 = t257 * t151 - t222;
t166 = t170 * qJD(3);
t94 = t142 * t122;
t85 = t94 * qJD(1);
t42 = t85 * pkin(2) + (-pkin(5) * t166 + (pkin(1) * t148 - t170 * pkin(5)) * qJD(2)) * qJD(1);
t245 = pkin(1) * qJD(2);
t141 = t147 * t245;
t124 = pkin(5) * t142 + t141;
t246 = pkin(1) * qJD(1);
t205 = t151 * t246;
t77 = pkin(2) * t114 - pkin(5) * t116 - t205;
t65 = t124 * t150 + t146 * t77;
t21 = -qJD(4) * t65 - t146 * t171 + t150 * t42;
t16 = t85 * pkin(3) + t21;
t64 = -t124 * t146 + t150 * t77;
t20 = qJD(4) * t64 + t146 * t42 + t150 * t171;
t237 = t149 * t65;
t111 = qJD(4) + t114;
t45 = pkin(3) * t111 + t64;
t27 = t145 * t45 + t237;
t5 = -t27 * qJD(5) - t145 * t20 + t149 * t16;
t125 = -t142 * pkin(2) - t185;
t78 = t97 * pkin(3) + t125;
t263 = t175 * t78 + t5;
t109 = qJD(5) + t111;
t155 = (t170 * qJD(2) + t166) * qJD(1);
t214 = qJD(4) * t146;
t201 = t116 * t213 + t142 * t214 + t146 * t155;
t212 = qJD(5) * t145;
t55 = t116 * t214 - t142 * t213 - t150 * t155;
t13 = t145 * t201 + t149 * t55 + t97 * t210 + t99 * t212;
t262 = t109 * t60 - t13;
t4 = (qJD(5) * t45 + t20) * t149 + t145 * t16 - t65 * t212;
t261 = t60 * t78 - t4;
t160 = t175 * qJD(5) + t145 * t55 - t149 * t201;
t260 = -t109 * t175 + t160;
t230 = t114 * t146;
t259 = (t214 + t230) * pkin(3);
t193 = t257 * qJD(3);
t93 = t142 * t222 + (-t194 - t193) * t151;
t256 = pkin(1) * t147;
t255 = t150 * pkin(3);
t253 = t99 * t97;
t156 = t208 * t258;
t110 = t116 * pkin(3);
t86 = pkin(2) * t116 + pkin(5) * t114;
t71 = -t146 * t185 + t150 * t86;
t56 = t110 + t71;
t72 = t146 * t86 + t150 * t185;
t252 = t156 * pkin(5) - t145 * t72 + t149 * t56;
t251 = pkin(5) * t265 + t145 * t56 + t149 * t72;
t138 = pkin(5) + t256;
t198 = t150 * t257;
t199 = t146 * t257;
t244 = pkin(1) * qJD(3);
t76 = pkin(1) * t218 + t86;
t58 = t150 * t76 + t110;
t250 = t149 * t58 - t76 * t225 - (-t145 * t198 - t149 * t199) * t244 + t156 * t138;
t249 = t145 * t58 + t76 * t223 - (-t145 * t199 + t149 * t198) * t244 + t265 * t138;
t243 = t111 * t97;
t242 = t111 * t99;
t241 = t145 * t65;
t240 = t146 * t55;
t239 = t146 * t85;
t238 = t146 * t93;
t236 = t150 * t64;
t235 = t150 * t85;
t234 = t150 * t99;
t19 = t20 * t150;
t69 = t85 * t170;
t233 = t97 * t146;
t232 = t109 * t116;
t231 = t111 * t116;
t187 = t111 * t150;
t229 = t116 * t114;
t228 = t122 * t150;
t226 = t125 * t114;
t221 = t148 * t151;
t219 = t148 ^ 2 - t151 ^ 2;
t216 = qJD(2) * t148;
t215 = qJD(3) * t147;
t211 = qJD(5) * t146;
t209 = qJD(1) * qJD(2);
t207 = -t114 * t236 - t65 * t230 + t19;
t204 = pkin(1) * t216;
t203 = pkin(1) * t215;
t153 = qJD(1) ^ 2;
t200 = t153 * t221;
t197 = t257 * t142;
t196 = t122 * t214;
t189 = -t64 * t116 + t125 * t214;
t188 = qJD(4) * t138 + t76;
t186 = pkin(1) * t193;
t184 = qJD(3) * t141;
t139 = -t257 * pkin(1) - pkin(2);
t183 = t209 * t221;
t182 = t203 + t259;
t181 = t201 * t150;
t180 = -t141 + t259;
t54 = pkin(2) * t94 + pkin(5) * t93 + t204;
t89 = -pkin(1) * t151 - pkin(2) * t170 - pkin(5) * t122;
t179 = pkin(3) * t94 + t150 * t54 + (-t211 - t214) * t89;
t178 = t65 * t116 + t125 * t213 + t146 * t184;
t176 = -0.2e1 * t183;
t174 = t146 * t65 + t236;
t173 = -t138 * t85 + t226;
t26 = t149 * t45 - t241;
t46 = t201 * pkin(3) + t184;
t169 = -t116 * t26 - t258 * t46 + t267 * t78;
t168 = t116 * t27 + t46 * t121 + t266 * t78;
t167 = -t5 * t121 + t258 * t4 - t266 * t26 - t267 * t27;
t165 = qJD(4) * t122 * t125 + t111 * t54 + t85 * t89;
t68 = -pkin(3) * t170 + t150 * t89;
t164 = qJD(5) * t68 + t146 * t54 + t89 * t213;
t161 = -t174 * qJD(4) - t21 * t146;
t159 = -qJD(4) * t111 * t89 + t122 * t184 - t125 * t93;
t158 = t161 + t19;
t154 = pkin(1) ^ 2;
t152 = qJD(2) ^ 2;
t140 = -pkin(2) - t255;
t126 = t139 - t255;
t118 = t258 * pkin(5);
t117 = t121 * pkin(5);
t105 = t116 * t205;
t104 = t114 * t205;
t102 = t258 * t138;
t101 = t121 * t138;
t83 = t258 * t122;
t82 = t121 * t122;
t70 = -t114 ^ 2 + t116 ^ 2;
t67 = t116 * t142 - t85;
t66 = t114 * t142 + t155;
t38 = t145 * t68 + t89 * t223;
t37 = t149 * t68 - t89 * t225;
t32 = t149 * t64 - t241;
t31 = -t145 * t64 - t237;
t30 = t111 * t187 - t116 * t99 + t239;
t29 = -t111 ^ 2 * t146 + t116 * t97 + t235;
t25 = t111 * t233 - t181;
t24 = t99 * t187 - t240;
t23 = -t93 * t224 + (t208 * t228 - t238) * t149 + (-t122 * t211 - t196) * t145;
t22 = t122 * t265 + t258 * t93;
t12 = -t109 * t267 + t116 * t60 + t258 * t85;
t11 = t109 * t266 + t116 * t175 + t121 * t85;
t10 = (-t55 - t243) * t150 + (-t201 - t242) * t146;
t9 = -t164 * t145 + t179 * t149;
t8 = t179 * t145 + t164 * t149;
t7 = t160 * t258 + t267 * t60;
t6 = -t121 * t13 - t175 * t266;
t1 = t121 * t160 - t13 * t258 + t175 * t267 - t266 * t60;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t183, -0.2e1 * t219 * t209, t152 * t151, t176, -t152 * t148, 0, 0, 0, 0, 0, -t116 * t93 + t122 * t155, t93 * t114 - t116 * t94 - t122 * t85 + t155 * t170, -t93 * t142, t114 * t94 - t69, -t94 * t142, 0, (t114 * t216 - t151 * t85 + (-t151 * t94 - t170 * t216) * qJD(1)) * pkin(1), t116 * t204 + (t122 * t216 + 0.2e1 * t151 * t93) * t246, (t257 * t93 - t147 * t94 + (t122 * t147 + t170 * t257) * qJD(3)) * t245, t154 * t176, -t93 * t234 + (-t55 * t150 - t99 * t214) * t122, (t99 * t146 + t150 * t97) * t93 + (-t181 + t240 + (t233 - t234) * qJD(4)) * t122, t85 * t228 + t170 * t55 + t94 * t99 + (-t150 * t93 - t196) * t111, -t93 * t233 + (t146 * t201 + t97 * t213) * t122, -t122 * t239 + t201 * t170 - t97 * t94 + (-t122 * t213 + t238) * t111, t111 * t94 - t69, t146 * t159 + t150 * t165 - t170 * t21 + t64 * t94, -t146 * t165 + t150 * t159 + t170 * t20 - t65 * t94, (-t21 * t122 - t54 * t99 + t89 * t55 + t64 * t93 + (-t122 * t65 - t89 * t97) * qJD(4)) * t150 + (-t54 * t97 - t89 * t201 - t20 * t122 + t65 * t93 + (t64 * t122 + t89 * t99) * qJD(4)) * t146, t174 * t54 + (t20 * t146 + t21 * t150 + (-t146 * t64 + t150 * t65) * qJD(4)) * t89, -t13 * t83 + t175 * t22, t13 * t82 + t160 * t83 + t175 * t23 + t22 * t60, -t109 * t22 + t13 * t170 - t175 * t94 + t83 * t85, -t160 * t82 + t23 * t60, -t109 * t23 - t160 * t170 - t60 * t94 - t82 * t85, t109 * t94 - t69, t109 * t9 - t170 * t5 + t23 * t78 + t26 * t94 + t37 * t85 + t46 * t82 + (-t60 * t238 + (-t146 * t160 + t213 * t60) * t122) * pkin(3), -t109 * t8 + t170 * t4 - t22 * t78 - t27 * t94 - t38 * t85 + t46 * t83 + (t175 * t238 + (-t13 * t146 - t175 * t213) * t122) * pkin(3), t13 * t37 + t160 * t38 + t175 * t9 + t22 * t26 - t23 * t27 - t4 * t82 - t5 * t83 - t60 * t8, t26 * t9 + t27 * t8 + t37 * t5 + t38 * t4 + (-t78 * t238 + (t146 * t46 + t213 * t78) * t122) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t219 * t153, 0, t200, 0, 0, 0, 0, 0, 0, t229, t70, t66, -t229, t67, 0, t105 + (-t114 * t218 + (-qJD(2) - t142) * t215) * pkin(1), -t104 + (-t116 * t218 + (-t194 - t197) * qJD(3)) * pkin(1), -t85 * t256 - pkin(1) * t142 * t195 * t170 + (t203 + t141) * t116 + (-t186 - t185) * t114, t154 * t200, t24, t10, t30, t25, t29, -t231, t139 * t201 + t173 * t146 - t188 * t187 + (-t111 * t199 + (-qJD(2) * t150 + t97) * t147) * t244 + t189, t99 * t203 - t139 * t55 + t173 * t150 + (t146 * t188 - t150 * t186) * t111 + t178, (-t97 * t186 - t138 * t201 + t76 * t99 + (t138 * t99 - t64) * qJD(4)) * t150 + (t99 * t186 - t138 * t55 + t76 * t97 - t21 + (t138 * t97 - t65) * qJD(4)) * t146 + t207, -t174 * t76 + t158 * t138 + (-t64 * t199 + t65 * t198 + (qJD(2) * t139 + t125) * t147) * t244, t6, t1, t11, t7, t12, -t232, -t101 * t85 - t250 * t109 - t126 * t160 + t182 * t60 + t169, -t102 * t85 + t249 * t109 - t126 * t13 - t175 * t182 + t168, -t101 * t13 + t102 * t160 - t175 * t250 + t249 * t60 + t167, -t101 * t5 + t102 * t4 + t126 * t46 + t182 * t78 - t249 * t27 - t250 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t70, t66, -t229, t67, 0, t105 + (-qJD(3) + t142) * t141, -t104 + (-t193 + t197) * t245, 0, 0, t24, t10, t30, t25, t29, -t231, -pkin(2) * t201 - t71 * t111 + t146 * t226 + (-t111 * t213 - t239) * pkin(5) + (-qJD(3) * t150 - t97) * t141 + t189, -t99 * t141 + t150 * t226 + pkin(2) * t55 + t111 * t72 + (t111 * t214 - t235) * pkin(5) + t178, t71 * t99 + t72 * t97 + (-t181 - t240 + (t233 + t234) * qJD(4)) * pkin(5) + t161 + t207, -t64 * t71 - t65 * t72 + (-pkin(2) * qJD(3) - t125) * t141 + t158 * pkin(5), t6, t1, t11, t7, t12, -t232, -t252 * t109 - t117 * t85 - t140 * t160 + t180 * t60 + t169, t251 * t109 - t118 * t85 - t13 * t140 - t175 * t180 + t168, -t117 * t13 + t118 * t160 - t175 * t252 + t251 * t60 + t167, -t117 * t5 + t118 * t4 + t140 * t46 + t180 * t78 - t251 * t27 - t252 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t97 ^ 2 + t99 ^ 2, -t55 + t243, -t253, -t201 + t242, t85, t65 * t111 - t125 * t99 + t21, t64 * t111 + t125 * t97 - t20, 0, 0, -t254, t264, t262, t254, t260, t85, -t109 * t31 + (-t109 * t212 + t149 * t85 - t60 * t99) * pkin(3) + t263, t109 * t32 + (-t109 * t210 - t145 * t85 + t175 * t99) * pkin(3) + t261, -t26 * t60 - t27 * t175 - t31 * t175 + t32 * t60 + (t13 * t149 + t160 * t145 + (-t145 * t175 - t149 * t60) * qJD(5)) * pkin(3), -t26 * t31 - t27 * t32 + (t145 * t4 + t149 * t5 - t78 * t99 + (-t145 * t26 + t149 * t27) * qJD(5)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, t264, t262, t254, t260, t85, t109 * t27 + t263, t109 * t26 + t261, 0, 0;];
tauc_reg  = t2;
