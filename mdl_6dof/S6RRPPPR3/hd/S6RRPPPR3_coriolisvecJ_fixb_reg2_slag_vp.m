% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:43
% EndTime: 2019-03-09 08:15:50
% DurationCPUTime: 3.28s
% Computational Cost: add. (4249->382), mult. (9240->508), div. (0->0), fcn. (5402->6), ass. (0->215)
t162 = sin(qJ(2));
t221 = t162 * qJD(1);
t129 = qJD(6) + t221;
t157 = cos(pkin(9));
t161 = sin(qJ(6));
t163 = cos(qJ(6));
t236 = t163 * t157;
t196 = t221 * t236;
t156 = sin(pkin(9));
t212 = t156 * t221;
t224 = qJD(6) * t163;
t237 = t161 * t156;
t255 = -qJD(6) * t237 + t157 * t224 - t161 * t212 + t196;
t273 = t255 * t129;
t164 = cos(qJ(2));
t269 = t236 - t237;
t82 = t269 * t164;
t227 = qJD(1) * t164;
t210 = t156 * t227;
t222 = t157 * qJD(2);
t97 = t210 - t222;
t209 = t157 * t227;
t223 = t156 * qJD(2);
t98 = t209 + t223;
t51 = t161 * t98 + t163 * t97;
t272 = t51 ^ 2;
t52 = t161 * t97 - t163 * t98;
t271 = t52 ^ 2;
t270 = t51 * t129;
t100 = t163 * t156 + t161 * t157;
t72 = t100 * t221;
t91 = t100 * qJD(6);
t256 = t72 + t91;
t208 = t256 * t129;
t137 = pkin(7) * t227;
t107 = -qJ(4) * t227 + t137;
t153 = qJD(2) * qJ(3);
t93 = -t107 - t153;
t154 = t162 ^ 2;
t155 = t164 ^ 2;
t268 = qJD(1) * (0.2e1 * t154 - t155);
t165 = -pkin(2) - pkin(3);
t213 = qJD(2) * t165;
t267 = -qJD(6) + t129;
t26 = qJD(2) * t72 + t52 * qJD(6);
t266 = t100 * t26 - t255 * t51;
t217 = qJD(1) * qJD(2);
t207 = t162 * t217;
t194 = t156 * t207;
t25 = t161 * (-qJD(6) * t98 + t194) - qJD(2) * t196 - t97 * t224;
t265 = t100 * t25 - t255 * t52;
t264 = t256 * t51 + t26 * t269;
t125 = qJ(4) * t207;
t152 = qJD(2) * qJD(3);
t219 = t164 * qJD(4);
t226 = qJD(2) * t162;
t178 = pkin(7) * t226 + t219;
t65 = qJD(1) * t178 - t125 - t152;
t149 = -qJ(5) + t165;
t176 = pkin(4) * t164 + t149 * t162;
t169 = qJD(2) * t176 + t164 * qJD(5);
t142 = t162 * qJD(3);
t206 = t164 * t217;
t230 = qJ(3) * t206 + qJD(1) * t142;
t35 = qJD(1) * t169 + t230;
t128 = pkin(7) * t206;
t225 = qJD(2) * t164;
t211 = qJ(4) * t225;
t220 = t162 * qJD(4);
t74 = t128 + (-t211 - t220) * qJD(1);
t66 = -qJD(2) * qJD(5) + t74;
t15 = -t156 * t66 + t157 * t35;
t239 = t157 * t162;
t182 = pkin(5) * t164 - pkin(8) * t239;
t175 = t182 * qJD(2);
t12 = qJD(1) * t175 + t15;
t16 = t156 * t35 + t157 * t66;
t13 = -pkin(8) * t194 + t16;
t188 = t161 * t12 + t163 * t13;
t189 = t162 * pkin(4) + t164 * qJ(5);
t95 = -qJD(1) * pkin(1) - pkin(2) * t227 - qJ(3) * t221;
t71 = pkin(3) * t227 + qJD(4) - t95;
t53 = t189 * qJD(1) + t71;
t132 = qJ(4) * t221;
t216 = pkin(7) * t221;
t190 = qJD(3) + t216;
t181 = -t132 + t190;
t69 = t149 * qJD(2) + t181;
t20 = -t156 * t69 + t157 * t53;
t14 = pkin(5) * t221 + t98 * pkin(8) + t20;
t21 = t156 * t53 + t157 * t69;
t17 = t97 * pkin(8) + t21;
t5 = t163 * t14 - t161 * t17;
t1 = t5 * qJD(6) + t188;
t6 = t161 * t14 + t163 * t17;
t2 = -t6 * qJD(6) + t163 * t12 - t161 * t13;
t263 = t1 * t269 - t2 * t100 - t255 * t5 - t256 * t6;
t134 = qJ(3) * t227;
t58 = qJD(1) * t176 + t134;
t33 = -t156 * t107 + t157 * t58;
t22 = qJD(1) * t182 + t33;
t34 = t157 * t107 + t156 * t58;
t27 = -pkin(8) * t212 + t34;
t257 = pkin(8) - t149;
t103 = t257 * t156;
t104 = t257 * t157;
t56 = t163 * t103 + t161 * t104;
t262 = -qJD(5) * t269 + qJD(6) * t56 - t161 * t22 - t163 * t27;
t57 = t161 * t103 - t163 * t104;
t261 = qJD(5) * t100 - qJD(6) * t57 + t161 * t27 - t163 * t22;
t259 = t52 * t51;
t258 = pkin(7) - qJ(4);
t160 = qJ(3) + pkin(4);
t229 = qJ(3) * t225 + t142;
t43 = t169 + t229;
t86 = t258 * t225 - t220;
t24 = t156 * t43 + t157 * t86;
t254 = qJD(2) * pkin(2);
t251 = t156 * t97;
t250 = t156 * t98;
t249 = t157 * t97;
t248 = t157 * t98;
t247 = t164 * t20;
t246 = t164 * t21;
t145 = t164 * qJ(4);
t117 = t164 * pkin(7) - t145;
t245 = t65 * t117;
t116 = t258 * t162;
t111 = -t164 * pkin(2) - t162 * qJ(3) - pkin(1);
t96 = t164 * pkin(3) - t111;
t68 = t189 + t96;
t42 = t157 * t116 + t156 * t68;
t243 = qJD(2) * t269;
t167 = qJD(1) ^ 2;
t141 = t154 * t167;
t242 = t155 * t167;
t241 = t156 * t162;
t240 = t156 * t164;
t238 = t157 * t164;
t235 = t164 * t167;
t166 = qJD(2) ^ 2;
t234 = t166 * t162;
t143 = t166 * t164;
t214 = pkin(5) * t156 - pkin(7);
t191 = t214 * t162;
t233 = -qJD(1) * t191 + qJD(3) - t132;
t232 = -qJD(4) - t71;
t105 = -t132 + t216;
t218 = qJD(3) + t105;
t215 = t157 * t141;
t23 = -t156 * t86 + t157 * t43;
t110 = t190 - t254;
t205 = -t110 - t254;
t197 = t162 * t213;
t60 = qJD(1) * t197 + t230;
t70 = t197 + t229;
t204 = qJD(1) * t70 + t60;
t203 = qJD(1) * t96 + t71;
t41 = -t156 * t116 + t157 * t68;
t202 = t232 * t162;
t201 = -0.2e1 * pkin(1) * t217;
t200 = qJD(1) * t111 + t95;
t199 = -t98 + t223;
t198 = t97 + t222;
t195 = -t100 * t206 - t273;
t193 = t162 * t206;
t187 = t15 * t157 + t16 * t156;
t186 = -t15 * t156 + t16 * t157;
t185 = -t156 * t21 - t157 * t20;
t184 = t156 * t20 - t157 * t21;
t30 = t162 * pkin(5) + pkin(8) * t238 + t41;
t36 = pkin(8) * t240 + t42;
t9 = -t161 * t36 + t163 * t30;
t10 = t161 * t30 + t163 * t36;
t183 = qJD(1) * t199;
t75 = pkin(2) * t207 - t230;
t87 = pkin(2) * t226 - t229;
t180 = -pkin(7) * t166 - qJD(1) * t87 - t75;
t179 = t25 * t269 + t256 * t52;
t177 = t198 * t221;
t173 = t100 * t226;
t78 = qJD(2) * pkin(4) + qJD(5) - t93;
t171 = -t149 * t225 + (qJD(2) * t160 + qJD(5) - t78) * t162;
t47 = pkin(5) * t194 - t65;
t109 = -pkin(7) * t207 + t152;
t115 = t137 + t153;
t168 = t109 * t164 + (t110 * t164 + (-t115 + t137) * t162) * qJD(2);
t151 = t157 ^ 2;
t150 = t156 ^ 2;
t146 = 0.2e1 * t152;
t135 = qJ(4) * t226;
t127 = t162 * t235;
t124 = t157 * pkin(5) + t160;
t123 = -t141 - t166;
t121 = t156 * t141;
t119 = -0.2e1 * t193;
t118 = 0.2e1 * t193;
t113 = -t141 + t242;
t106 = pkin(2) * t221 - t134;
t102 = (-t154 + t155) * t217;
t94 = 0.2e1 * t102;
t89 = -t214 * t164 - t145;
t85 = -t135 + t178;
t84 = t165 * t221 + t134;
t81 = t100 * t164;
t79 = t213 + t181;
t64 = qJD(2) * t191 + t135 - t219;
t46 = -t97 * pkin(5) + t78;
t40 = qJD(6) * t82 - t173;
t39 = t164 * t91 + t226 * t269;
t19 = -t162 * pkin(8) * t223 + t24;
t18 = t175 + t23;
t4 = -t10 * qJD(6) - t161 * t19 + t163 * t18;
t3 = t9 * qJD(6) + t161 * t18 + t163 * t19;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t94, t143, t119, -t234, 0, -pkin(7) * t143 + t162 * t201, pkin(7) * t234 + t164 * t201, 0, 0, t118, t143, -0.2e1 * t102, 0, t234, t119, t180 * t164 + t200 * t226, t168, t180 * t162 - t200 * t225, pkin(7) * t168 + t75 * t111 + t95 * t87, t119, t94, -t234, t118, t143, 0, t204 * t162 + (t203 * t164 - t85) * qJD(2), -t204 * t164 + (t162 * t203 + t86) * qJD(2), -t74 * t162 + t65 * t164 + (-t162 * t93 - t164 * t79) * qJD(2) + (-t162 * t86 + t164 * t85 + (-t116 * t164 + t117 * t162) * qJD(2)) * qJD(1), t74 * t116 + t60 * t96 + t71 * t70 + t79 * t86 + t93 * t85 - t245 (-t151 * t227 - t248) * t226 (t249 + (t98 + 0.2e1 * t209) * t156) * t226 (t157 * t268 - t164 * t98) * qJD(2) (-t150 * t227 - t251) * t226 (-t156 * t268 + t164 * t97) * qJD(2), t118, t65 * t240 + t85 * t97 + (qJD(1) * t23 + t15) * t162 + (t78 * t241 + t247 + (t117 * t241 + t164 * t41) * qJD(1)) * qJD(2), t65 * t238 + t85 * t98 + (-qJD(1) * t24 - t16) * t162 + (t78 * t239 - t246 + (t117 * t239 - t164 * t42) * qJD(1)) * qJD(2), t23 * t98 + t24 * t97 + t187 * t164 + ((-t156 * t42 - t157 * t41) * qJD(1) + t185) * t226, t15 * t41 + t16 * t42 + t20 * t23 + t21 * t24 - t78 * t85 - t245, t25 * t82 + t52 * t39, -t25 * t81 + t82 * t26 + t39 * t51 + t52 * t40, t39 * t129 - t25 * t162 + (-qJD(1) * t82 + t52) * t225, -t26 * t81 + t40 * t51, t40 * t129 - t26 * t162 + (qJD(1) * t81 + t51) * t225 (t129 + t221) * t225, t4 * t129 + t2 * t162 + t89 * t26 - t46 * t40 - t47 * t81 - t64 * t51 + (qJD(1) * t9 + t5) * t225, -t1 * t162 - t3 * t129 - t89 * t25 + t46 * t39 - t47 * t82 + t64 * t52 + (-qJD(1) * t10 - t6) * t225, t1 * t81 - t10 * t26 + t2 * t82 + t9 * t25 + t3 * t51 - t5 * t39 - t4 * t52 + t6 * t40, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + t46 * t64 + t47 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t113, 0, t127, 0, 0, t167 * pkin(1) * t162, pkin(1) * t235, 0, 0, -t127, 0, t113, 0, 0, t127 (t106 * t164 - t162 * t95) * qJD(1) ((t115 - t153) * t162 + (qJD(3) + t205) * t164) * qJD(1), t146 + (t106 * t162 + t164 * t95) * qJD(1), t109 * qJ(3) + t115 * qJD(3) - t95 * t106 + (t115 * t162 + t205 * t164) * qJD(1) * pkin(7), t127, -t113, 0, -t127, 0, 0, t105 * qJD(2) + t125 + t146 + (t232 * t164 + (-pkin(7) * qJD(2) - t84) * t162) * qJD(1), -t107 * qJD(2) + t128 + ((-qJ(4) * qJD(2) + t84) * t164 + t202) * qJD(1) (-t218 + t79 - t213) * t227, -t65 * qJ(3) - t79 * t107 + t74 * t165 - t218 * t93 - t71 * t84, -t183 * t239 (-t250 - t249 + (t150 - t151) * qJD(2)) * t221, -t164 * t183 - t215, t156 * t177, -t198 * t227 + t121, -t127, -t65 * t157 - t218 * t97 + (t156 * t171 - t162 * t33 - t247) * qJD(1), t65 * t156 - t218 * t98 + (t157 * t171 + t162 * t34 + t246) * qJD(1), -t33 * t98 - t34 * t97 + (-qJD(5) * t97 + t20 * t221 - t16) * t157 + (qJD(5) * t98 + t21 * t221 + t15) * t156, qJD(5) * t184 + t149 * t186 - t65 * t160 - t20 * t33 - t21 * t34 + t218 * t78, t265, t179 + t266, -t227 * t52 + t195, t264, t208 + (-t51 - t243) * t227, -t129 * t227, t124 * t26 + t47 * t269 - t233 * t51 - t256 * t46 + t261 * t129 + (qJD(2) * t56 - t5) * t227, -t47 * t100 - t124 * t25 + t233 * t52 - t255 * t46 - t262 * t129 + (-qJD(2) * t57 + t6) * t227, t56 * t25 - t57 * t26 - t261 * t52 + t262 * t51 - t263, t1 * t57 + t47 * t124 + t2 * t56 + t233 * t46 + t261 * t5 + t262 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, 0, t123, -t115 * qJD(2) + t95 * t221 + t128, 0, 0, 0, 0, 0, 0, t123, t127, 0, t93 * qJD(2) + t128 + (t202 - t211) * qJD(1), 0, 0, 0, 0, 0, 0, -t215 + (t97 - t210) * qJD(2), t121 + (t98 - t209) * qJD(2) (-t248 - t251) * t221, -t78 * qJD(2) + t185 * t221 + t186, 0, 0, 0, 0, 0, 0, qJD(2) * t51 + t195, t208 + (-t227 * t269 - t52) * qJD(2), -t264 - t265, -t46 * qJD(2) + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t206, 0.2e1 * t207, -t141 - t242 (-t164 * t93 + (t79 + t213) * t162) * qJD(1) + t230, 0, 0, 0, 0, 0, 0, -t121 + (-t97 + t222) * t227, -t215 + (-t98 - t223) * t227 (-t250 + t249 + (-t150 - t151) * qJD(2)) * t221 (-t162 * t184 + t164 * t78) * qJD(1) + t187, 0, 0, 0, 0, 0, 0, -t208 + (-t51 + t243) * t227, -t273 + (-qJD(2) * t100 + t52) * t227, t179 - t266, t1 * t100 + t2 * t269 + t227 * t46 + t255 * t6 - t256 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199 * t221, t177, -t97 ^ 2 - t98 ^ 2, -t20 * t98 - t21 * t97 - t65, 0, 0, 0, 0, 0, 0, t52 * t129 + t26, -t25 + t270, -t271 - t272, t5 * t52 - t6 * t51 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t271 - t272, -t25 - t270, t259, -qJD(1) * t173 + t267 * t52, t206, t6 * t129 - t46 * t52 + t2, t267 * t5 - t46 * t51 - t188, 0, 0;];
tauc_reg  = t7;
