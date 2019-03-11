% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:12
% EndTime: 2019-03-08 21:54:18
% DurationCPUTime: 3.03s
% Computational Cost: add. (4299->302), mult. (11310->444), div. (0->0), fcn. (9244->12), ass. (0->191)
t164 = cos(qJ(6));
t214 = qJD(6) * t164;
t156 = sin(pkin(12));
t158 = cos(pkin(12));
t162 = sin(qJ(3));
t166 = cos(qJ(3));
t135 = -t156 * t162 + t158 * t166;
t128 = t135 * qJD(2);
t165 = cos(qJ(5));
t116 = t165 * t128;
t136 = t156 * t166 + t158 * t162;
t130 = t136 * qJD(2);
t161 = sin(qJ(5));
t83 = -t161 * t130 + t116;
t274 = t164 * t83;
t281 = t214 - t274;
t153 = qJD(3) + qJD(5);
t234 = t83 * t153;
t129 = t136 * qJD(3);
t119 = qJD(2) * t129;
t213 = qJD(2) * qJD(3);
t201 = t166 * t213;
t202 = t162 * t213;
t120 = -t156 * t202 + t158 * t201;
t216 = qJD(5) * t161;
t43 = qJD(5) * t116 - t161 * t119 + t165 * t120 - t130 * t216;
t280 = t43 - t234;
t245 = -qJ(4) - pkin(8);
t197 = qJD(3) * t245;
t124 = t166 * qJD(4) + t162 * t197;
t125 = -t162 * qJD(4) + t166 * t197;
t167 = cos(qJ(2));
t157 = sin(pkin(6));
t221 = qJD(1) * t157;
t204 = t167 * t221;
t231 = -t156 * t124 + t158 * t125 + t136 * t204;
t230 = t158 * t124 + t156 * t125 - t135 * t204;
t223 = -qJD(6) + t83;
t279 = qJD(6) + t223;
t160 = sin(qJ(6));
t183 = t161 * t128 + t165 * t130;
t215 = qJD(6) * t160;
t24 = t153 * t214 + t164 * t43 - t183 * t215;
t71 = t160 * t153 + t164 * t183;
t25 = qJD(6) * t71 + t160 * t43;
t69 = -t164 * t153 + t160 * t183;
t278 = -t160 * t25 + t24 * t164 - t281 * t69;
t22 = t24 * t160;
t277 = t281 * t71 + t22;
t44 = qJD(5) * t183 + t165 * t119 + t161 * t120;
t37 = t160 * t44;
t72 = t223 * t214;
t243 = t37 - t72;
t248 = t71 * t183;
t276 = t223 * t274 + t243 - t248;
t251 = t130 * pkin(9);
t163 = sin(qJ(2));
t205 = t163 * t221;
t140 = qJD(2) * pkin(8) + t205;
t192 = qJ(4) * qJD(2) + t140;
t159 = cos(pkin(6));
t220 = qJD(1) * t159;
t203 = t162 * t220;
t99 = t192 * t166 + t203;
t91 = t156 * t99;
t147 = t166 * t220;
t98 = -t192 * t162 + t147;
t95 = qJD(3) * pkin(3) + t98;
t52 = t158 * t95 - t91;
t40 = qJD(3) * pkin(4) - t251 + t52;
t252 = t128 * pkin(9);
t241 = t158 * t99;
t53 = t156 * t95 + t241;
t45 = t53 + t252;
t16 = -t161 * t45 + t165 * t40;
t14 = -t153 * pkin(5) - t16;
t275 = t14 * t83;
t273 = t183 * t83;
t132 = t135 * qJD(3);
t272 = t132 * pkin(9) - t231;
t271 = t129 * pkin(9) - t230;
t270 = t160 * t223;
t235 = t183 * t153;
t269 = -t44 + t235;
t217 = qJD(3) * t162;
t268 = pkin(3) * t217 - t205;
t266 = t183 ^ 2 - t83 ^ 2;
t49 = pkin(5) * t183 - t83 * pkin(10);
t185 = qJD(4) + t204;
t67 = (-t162 * t140 + t147) * qJD(3) + (-qJ(4) * t217 + t185 * t166) * qJD(2);
t68 = (-t166 * t140 - t203) * qJD(3) + (-qJD(3) * t166 * qJ(4) - t185 * t162) * qJD(2);
t30 = -t156 * t67 + t158 * t68;
t27 = -t120 * pkin(9) + t30;
t31 = t156 * t68 + t158 * t67;
t28 = -t119 * pkin(9) + t31;
t2 = (qJD(5) * t40 + t28) * t165 + t161 * t27 - t45 * t216;
t208 = -t166 * pkin(3) - pkin(2);
t121 = t208 * qJD(2) + qJD(4) - t204;
t86 = -t128 * pkin(4) + t121;
t265 = -t86 * t83 - t2;
t247 = t183 * t69;
t262 = t223 * t183;
t39 = t164 * t44;
t261 = -t215 * t223 - t39;
t187 = t129 * pkin(4) + t268;
t17 = t161 * t40 + t165 * t45;
t15 = t153 * pkin(10) + t17;
t29 = -pkin(5) * t83 - pkin(10) * t183 + t86;
t184 = t160 * t15 - t164 * t29;
t260 = t14 * t215 + t183 * t184;
t3 = t17 * qJD(5) + t161 * t28 - t165 * t27;
t5 = t164 * t15 + t160 * t29;
t259 = t14 * t214 + t3 * t160 + t5 * t183;
t182 = t165 * t135 - t161 * t136;
t142 = t245 * t162;
t143 = t245 * t166;
t102 = t158 * t142 + t156 * t143;
t75 = -t136 * pkin(9) + t102;
t103 = t156 * t142 - t158 * t143;
t76 = t135 * pkin(9) + t103;
t35 = t161 * t76 - t165 * t75;
t255 = t35 * qJD(5) + t272 * t161 + t271 * t165;
t110 = -t135 * pkin(4) + t208;
t90 = t161 * t135 + t165 * t136;
t33 = -pkin(5) * t182 - t90 * pkin(10) + t110;
t36 = t161 * t75 + t165 * t76;
t50 = qJD(5) * t182 - t161 * t129 + t165 * t132;
t258 = (qJD(6) * t29 + t2) * t182 + t14 * t50 + t3 * t90 - (-qJD(6) * t33 + t255) * t223 - t36 * t44;
t257 = -t183 * t86 - t3;
t254 = t36 * qJD(5) - t271 * t161 + t272 * t165;
t253 = pkin(3) * t156;
t250 = t14 * t90;
t249 = t33 * t44;
t246 = t90 * t44;
t57 = t158 * t98 - t91;
t242 = qJD(2) * pkin(2);
t239 = t160 * t71;
t150 = t158 * pkin(3) + pkin(4);
t179 = t165 * t150 - t161 * t253;
t55 = -t156 * t98 - t241;
t46 = t55 - t252;
t47 = t57 - t251;
t233 = -t179 * qJD(5) + t161 * t46 + t165 * t47;
t180 = t161 * t150 + t165 * t253;
t232 = t180 * qJD(5) - t161 * t47 + t165 * t46;
t228 = t157 * t163;
t227 = t157 * t167;
t169 = qJD(2) ^ 2;
t226 = t157 * t169;
t168 = qJD(3) ^ 2;
t225 = t168 * t162;
t224 = t168 * t166;
t126 = pkin(3) * t202 + qJD(2) * t205;
t222 = t162 ^ 2 - t166 ^ 2;
t219 = qJD(2) * t162;
t218 = qJD(2) * t163;
t210 = t90 * t215;
t209 = t163 * t226;
t207 = t157 * t218;
t206 = qJD(2) * t227;
t104 = pkin(3) * t219 + t130 * pkin(4);
t123 = pkin(10) + t180;
t194 = qJD(6) * t123 + t104 + t49;
t191 = t162 * t206;
t190 = t166 * t206;
t85 = t119 * pkin(4) + t126;
t51 = qJD(5) * t90 + t165 * t129 + t161 * t132;
t189 = t51 * pkin(5) - t50 * pkin(10) + t187;
t188 = -t223 * t50 + t246;
t133 = t159 * t166 - t162 * t228;
t134 = t159 * t162 + t166 * t228;
t78 = t158 * t133 - t156 * t134;
t79 = t156 * t133 + t158 * t134;
t41 = t161 * t79 - t165 * t78;
t42 = t161 * t78 + t165 * t79;
t181 = -t270 * t83 - t261;
t178 = -t160 * t42 - t164 * t227;
t177 = t160 * t227 - t164 * t42;
t176 = t242 * qJD(2);
t175 = -t123 * t44 - t223 * t233 - t275;
t174 = -0.2e1 * qJD(3) * t242;
t122 = -pkin(5) - t179;
t97 = -qJD(3) * t134 - t191;
t96 = qJD(3) * t133 + t190;
t56 = t156 * t97 + t158 * t96;
t54 = -t156 * t96 + t158 * t97;
t11 = t44 * pkin(5) - t43 * pkin(10) + t85;
t10 = t164 * t11;
t7 = t42 * qJD(5) + t161 * t56 - t165 * t54;
t6 = -t41 * qJD(5) + t161 * t54 + t165 * t56;
t1 = [0, 0, -t209, -t167 * t226, 0, 0, 0, 0, 0, -t166 * t209 + (t97 - t191) * qJD(3), t162 * t209 + (-t96 - t190) * qJD(3), -t79 * t119 - t78 * t120 + t56 * t128 - t54 * t130, t30 * t78 + t31 * t79 + t52 * t54 + t53 * t56 + (t121 * t218 - t126 * t167) * t157, 0, 0, 0, 0, 0, -t7 * t153 + (-t167 * t44 - t218 * t83) * t157, -t6 * t153 + (-t167 * t43 + t183 * t218) * t157, 0, 0, 0, 0, 0 -(qJD(6) * t177 - t160 * t6 + t164 * t207) * t223 + t178 * t44 + t7 * t69 + t41 * t25 (qJD(6) * t178 + t160 * t207 + t164 * t6) * t223 + t177 * t44 + t7 * t71 + t41 * t24; 0, 0, 0, 0, 0.2e1 * t162 * t201, -0.2e1 * t222 * t213, t224, -t225, 0, -pkin(8) * t224 + t162 * t174, pkin(8) * t225 + t166 * t174, -t102 * t120 - t103 * t119 + t230 * t128 - t53 * t129 - t231 * t130 - t52 * t132 + t31 * t135 - t30 * t136, t30 * t102 + t31 * t103 + t268 * t121 + t126 * t208 + t230 * t53 + t231 * t52, t183 * t50 + t43 * t90, t182 * t43 - t183 * t51 + t50 * t83 - t246, t50 * t153, -t51 * t153, 0, t110 * t44 - t254 * t153 - t182 * t85 - t187 * t83 + t86 * t51, t110 * t43 + t255 * t153 + t183 * t187 + t86 * t50 + t85 * t90, -t71 * t210 + (t24 * t90 + t50 * t71) * t164 (-t164 * t69 - t239) * t50 + (-t22 - t164 * t25 + (t160 * t69 - t164 * t71) * qJD(6)) * t90, t164 * t188 - t182 * t24 + t210 * t223 + t71 * t51, -t160 * t188 + t182 * t25 - t69 * t51 + t72 * t90, -t182 * t44 - t223 * t51, -t10 * t182 + t35 * t25 - t184 * t51 + t254 * t69 + (t249 - t189 * t223 + (t15 * t182 + t223 * t36 + t250) * qJD(6)) * t164 + t258 * t160, t35 * t24 - t5 * t51 + t254 * t71 + (-t249 + (-qJD(6) * t15 + t11) * t182 - qJD(6) * t250 - (qJD(6) * t36 - t189) * t223) * t160 + t258 * t164; 0, 0, 0, 0, -t162 * t169 * t166, t222 * t169, 0, 0, 0, t162 * t176, t166 * t176 (t53 + t55) * t130 + (t52 - t57) * t128 + (-t119 * t156 - t120 * t158) * pkin(3), -t52 * t55 - t53 * t57 + (-t121 * t219 + t156 * t31 + t158 * t30) * pkin(3), -t273, t266, t280, t269, 0, t104 * t83 - t153 * t232 + t257, -t104 * t183 + t153 * t233 + t265, t277, t223 * t239 + t278, t276, t181 + t247, t262, t122 * t25 + t232 * t69 + (t194 * t223 - t3) * t164 + t175 * t160 + t260, t122 * t24 + t164 * t175 - t194 * t270 + t232 * t71 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 ^ 2 - t130 ^ 2, -t53 * t128 + t52 * t130 + t126, 0, 0, 0, 0, 0, t44 + t235, t43 + t234, 0, 0, 0, 0, 0, t181 - t247, -t164 * t223 ^ 2 - t248 - t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t266, t280, t269, 0, t17 * t153 + t257, t16 * t153 + t265, t277, t270 * t71 + t278, t276, -t223 * t270 + t247 + t39, t262, -pkin(5) * t25 - t3 * t164 + (-t160 * t16 + t164 * t49) * t223 - t17 * t69 - t160 * t275 - t243 * pkin(10) + t260, -pkin(5) * t24 - (t164 * t16 + t160 * t49) * t223 - t17 * t71 - t14 * t274 + t261 * pkin(10) + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t69 ^ 2 + t71 ^ 2, -t223 * t69 + t24, -t223 * t71 - t25, t44, -t14 * t71 - t160 * t2 - t279 * t5 + t10, -t160 * t11 + t14 * t69 - t164 * t2 + t279 * t184;];
tauc_reg  = t1;
