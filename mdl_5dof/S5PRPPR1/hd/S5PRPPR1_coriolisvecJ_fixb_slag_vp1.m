% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:46
% EndTime: 2019-12-05 15:21:57
% DurationCPUTime: 5.04s
% Computational Cost: add. (8260->377), mult. (8092->536), div. (0->0), fcn. (7638->8), ass. (0->192)
t163 = pkin(9) + qJ(5);
t159 = sin(t163);
t161 = cos(t163);
t166 = sin(pkin(8));
t168 = cos(pkin(8));
t164 = pkin(7) + qJ(2);
t162 = cos(t164);
t160 = sin(t164);
t247 = t160 * t168;
t104 = t159 * t247 + t161 * t162;
t105 = -t162 * t159 + t161 * t247;
t248 = t160 * t166;
t57 = Icges(6,5) * t105 - Icges(6,6) * t104 + Icges(6,3) * t248;
t93 = Icges(6,4) * t105;
t61 = Icges(6,2) * t104 - Icges(6,6) * t248 - t93;
t92 = Icges(6,4) * t104;
t63 = Icges(6,1) * t105 + Icges(6,5) * t248 - t92;
t23 = (t159 * t61 + t161 * t63) * t166 - t168 * t57;
t289 = rSges(5,3) + qJ(4);
t180 = (-rSges(6,1) * t159 - rSges(6,2) * t161) * t166;
t114 = qJD(5) * t180;
t226 = qJD(5) * t166;
t288 = t114 * t226;
t287 = t104 * t61 + t105 * t63;
t244 = t168 * t162;
t106 = -t159 * t244 + t160 * t161;
t107 = t159 * t160 + t161 * t244;
t286 = -t106 * t61 + t107 * t63;
t245 = t162 * t166;
t283 = t57 * t245;
t165 = sin(pkin(9));
t167 = cos(pkin(9));
t246 = t162 * t165;
t196 = rSges(5,2) * (t162 * t167 + t165 * t247) - rSges(5,1) * (t167 * t247 - t246);
t249 = t160 * t165;
t282 = -(t167 * t244 + t249) * rSges(5,1) - (t160 * t167 - t165 * t244) * rSges(5,2);
t150 = -qJD(5) * t168 + qJD(2);
t103 = -rSges(6,3) * t168 + (rSges(6,1) * t161 - rSges(6,2) * t159) * t166;
t212 = t103 * t226;
t195 = rSges(6,1) * t105 - rSges(6,2) * t104;
t66 = rSges(6,3) * t248 + t195;
t281 = -t150 * t66 + t160 * t212;
t16 = t283 + t286;
t278 = t16 - t283;
t131 = t162 * pkin(2) + t160 * qJ(3);
t158 = pkin(4) * t167 + pkin(3);
t169 = -pkin(6) - qJ(4);
t277 = pkin(4) * t249 + t158 * t244 - t169 * t245 + t131;
t187 = rSges(4,1) * t244 - rSges(4,2) * t245 + t160 * rSges(4,3);
t276 = t131 + t187;
t275 = pkin(3) * t244 + t245 * t289 + t131 - t282;
t59 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t245;
t253 = Icges(6,4) * t107;
t62 = Icges(6,2) * t106 + Icges(6,6) * t245 + t253;
t94 = Icges(6,4) * t106;
t65 = Icges(6,1) * t107 + Icges(6,5) * t245 + t94;
t15 = -t104 * t62 + t105 * t65 + t59 * t248;
t274 = t160 * (-Icges(6,2) * t105 + t63 - t92) + t162 * (-Icges(6,2) * t107 + t65 + t94);
t230 = qJD(2) * t166;
t216 = t160 * t230;
t74 = qJD(2) * t104 - qJD(5) * t107;
t75 = -qJD(2) * t105 + qJD(5) * t106;
t266 = t75 * rSges(6,1) + t74 * rSges(6,2);
t39 = -rSges(6,3) * t216 + t266;
t76 = qJD(2) * t106 - qJD(5) * t105;
t77 = qJD(2) * t107 - qJD(5) * t104;
t199 = -rSges(6,1) * t77 - rSges(6,2) * t76;
t215 = t162 * t230;
t40 = rSges(6,3) * t215 - t199;
t192 = -t160 * t39 + t162 * t40;
t68 = t107 * rSges(6,1) + t106 * rSges(6,2) + rSges(6,3) * t245;
t9 = ((-t160 * t66 - t162 * t68) * qJD(2) + t192) * t226;
t273 = m(6) * t9;
t272 = t160 / 0.2e1;
t271 = -t162 / 0.2e1;
t270 = pkin(3) * t168;
t269 = pkin(4) * t165;
t265 = t196 * qJD(2);
t194 = qJ(4) * t166 + t270;
t228 = qJD(4) * t166;
t138 = t162 * t228;
t232 = qJD(2) * t160;
t225 = qJD(2) * qJD(3);
t151 = qJD(3) * t160;
t231 = qJD(2) * t162;
t234 = qJ(3) * t231 + t151;
t254 = t160 * t225 + qJD(2) * (-pkin(2) * t232 + t234);
t205 = t254 + (-t194 * t232 + 0.2e1 * t138) * qJD(2);
t221 = qJD(2) * t269;
t238 = t162 * t221 + t169 * t216;
t250 = t158 * t168;
t10 = -t162 * t288 + t150 * t39 + ((t212 + (t194 - t250) * qJD(2)) * t160 + t238) * qJD(2) + t205;
t262 = t10 * t162;
t147 = t162 * t225;
t209 = (pkin(3) - t158) * t168;
t227 = qJD(5) * t103;
t243 = qJ(4) + t169;
t152 = qJD(3) * t162;
t108 = qJD(2) * t131 - t152;
t214 = t160 * t228;
t255 = -t194 * t231 - t108 - t214;
t11 = t160 * t288 - t150 * t40 + t147 + ((-t221 - t228) * t160 + (qJD(2) * t209 + (qJD(2) * t243 + t227) * t166) * t162 + t255) * qJD(2);
t261 = t11 * t160;
t14 = t248 * t57 + t287;
t260 = t14 * t160;
t100 = -Icges(6,3) * t168 + (Icges(6,5) * t161 - Icges(6,6) * t159) * t166;
t251 = Icges(6,4) * t161;
t101 = -Icges(6,6) * t168 + (-Icges(6,2) * t159 + t251) * t166;
t252 = Icges(6,4) * t159;
t102 = -Icges(6,5) * t168 + (Icges(6,1) * t161 - t252) * t166;
t31 = t100 * t248 - t101 * t104 + t102 * t105;
t259 = t150 * t31;
t236 = t138 + t151;
t123 = t194 * t160;
t154 = t162 * qJ(3);
t129 = pkin(2) * t160 - t154;
t240 = -t123 - t129;
t224 = pkin(4) * t246;
t73 = t224 + (t166 * t243 + t209) * t160;
t21 = (t73 + t240) * qJD(2) + t236 + t281;
t258 = t162 * t21;
t256 = -rSges(6,3) + t169;
t122 = (-Icges(6,1) * t159 - t251) * t166;
t242 = -t101 + t122;
t121 = (-Icges(6,2) * t161 - t252) * t166;
t241 = t102 + t121;
t237 = rSges(4,2) * t216 + rSges(4,3) * t231;
t235 = rSges(4,2) * t248 + t162 * rSges(4,3);
t233 = -qJD(2) * t129 + t151;
t229 = qJD(4) * t160;
t17 = t106 * t62 + t107 * t65 + t59 * t245;
t223 = rSges(4,1) * t247;
t218 = t138 + t234;
t217 = -pkin(2) - t270;
t210 = -rSges(4,1) * t168 - pkin(2);
t208 = -t226 / 0.2e1;
t207 = t226 / 0.2e1;
t206 = -pkin(2) - t250;
t204 = -qJD(2) * t123 + t138 + t233;
t203 = t160 * t208;
t202 = t160 * t207;
t201 = t162 * t208;
t200 = t162 * t207;
t197 = t282 * qJD(2);
t22 = t150 * t68 - t152 + (-t162 * t227 + t229) * t166 + t277 * qJD(2);
t193 = t160 * t21 - t162 * t22;
t191 = -t160 * t68 + t162 * t66;
t190 = t160 * (-Icges(6,5) * t104 - Icges(6,6) * t105) + t162 * (Icges(6,5) * t106 - Icges(6,6) * t107);
t189 = qJD(2) * t203;
t188 = qJD(2) * t200;
t181 = t166 * t190;
t179 = (t15 * t162 + t260) * t166;
t178 = (t16 * t160 + t162 * t17) * t166;
t120 = (-Icges(6,5) * t159 - Icges(6,6) * t161) * t166;
t176 = (Icges(6,1) * t106 - t253 - t62) * t162 + (-Icges(6,1) * t104 + t61 - t93) * t160;
t174 = -rSges(5,3) * t248 + t196;
t33 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t216;
t34 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t215;
t35 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t216;
t36 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t215;
t37 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t216;
t38 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t215;
t173 = ((t106 * t36 + t107 * t38 - t61 * t74 + t63 * t75 + (t162 * t34 - t232 * t57) * t166) * t160 + t162 * (t106 * t35 + t107 * t37 + t62 * t74 + t65 * t75 + (t162 * t33 - t232 * t59) * t166) + (t16 * t162 - t17 * t160) * qJD(2)) * t166;
t172 = (t160 * (-t104 * t36 + t105 * t38 - t61 * t76 + t63 * t77 + (t160 * t34 + t231 * t57) * t166) + t162 * (-t104 * t35 + t105 * t37 + t62 * t76 + t65 * t77 + (t160 * t33 + t231 * t59) * t166) + (t14 * t162 - t15 * t160) * qJD(2)) * t166;
t24 = -t168 * t59 + (-t159 * t62 + t161 * t65) * t166;
t7 = -t168 * t34 + (-t159 * t36 + t161 * t38 + (-t159 * t63 + t161 * t61) * qJD(5)) * t166;
t8 = -t168 * t33 + (-t159 * t35 + t161 * t37 + (-t159 * t65 - t161 * t62) * qJD(5)) * t166;
t171 = (t160 * t7 + t162 * t8 + (-t160 * t24 + t23 * t162) * qJD(2)) * t166;
t113 = qJD(5) * t122;
t112 = qJD(5) * t121;
t111 = qJD(5) * t120;
t91 = t223 - t235;
t87 = qJD(2) * t276 - t152;
t86 = t151 + (-t129 - t91) * qJD(2);
t85 = rSges(6,1) * t106 - rSges(6,2) * t107;
t84 = -rSges(6,1) * t104 - rSges(6,2) * t105;
t56 = t147 + (-qJD(2) * t187 - t108) * qJD(2);
t55 = qJD(2) * (-qJD(2) * t223 + t237) + t254;
t45 = qJD(2) * t275 - t152 + t214;
t44 = (t174 + t240) * qJD(2) + t236;
t32 = t100 * t245 + t101 * t106 + t102 * t107;
t30 = t32 * t150;
t29 = t147 + ((-rSges(5,3) * t231 - t229) * t166 + t197 + t255) * qJD(2);
t28 = qJD(2) * (-rSges(5,3) * t216 + t265) + t205;
t27 = -qJD(4) * t168 + t191 * t226 + qJD(1);
t26 = -t111 * t168 + (-t112 * t159 + t113 * t161 + (-t101 * t161 - t102 * t159) * qJD(5)) * t166;
t25 = t26 * t150;
t13 = t101 * t76 + t102 * t77 - t104 * t112 + t105 * t113 + (t100 * t231 + t111 * t160) * t166;
t12 = t101 * t74 + t102 * t75 + t106 * t112 + t107 * t113 + (-t100 * t232 + t111 * t162) * t166;
t6 = qJD(5) * t178 + t30;
t5 = qJD(5) * t179 + t259;
t1 = [t273; t25 + (t30 + ((t14 + t17 - t287) * t162 + t278 * t160) * t226) * t203 + (t8 + t12) * t200 + (t24 + t32) * t189 + (t23 + t31) * t188 + (t11 * (t154 - t195 + t224) + t21 * (t152 + t199) + t10 * (t68 + t277) + t22 * (t218 + t238 + t266) + (t11 * t206 + (-t21 * qJD(4) + t11 * t256) * t166) * t160 + ((t166 * t256 + t206) * t258 + (t21 * (-qJ(3) - t269) + t22 * (-rSges(6,3) * t166 + t206)) * t160) * qJD(2) - (qJD(2) * t73 + t204 - t21 + t281) * t22) * m(6) + (t29 * (t154 + t196) + t44 * (t152 + t197) + t28 * t275 + t45 * (t218 + t265) + (t29 * t217 + (-t44 * qJD(4) - t289 * t29) * t166) * t160 + (t44 * (-t166 * t289 + t217) * t162 + (-t44 * qJ(3) + t45 * (-rSges(5,3) * t166 - pkin(2) - t194)) * t160) * qJD(2) - (qJD(2) * t174 + t204 - t44) * t45) * m(5) + (t56 * (t160 * t210 + t154 + t235) + t86 * t152 + t55 * t276 + t87 * (t234 + t237) + (t86 * (rSges(4,2) * t166 + t210) * t162 + (t86 * (-rSges(4,3) - qJ(3)) + t87 * t210) * t160) * qJD(2) - (-qJD(2) * t91 + t233 - t86) * t87) * m(4) + (t6 + t7 + t13) * t202 + (t5 - t259 + ((-t15 + t278 - t286) * t162 - t260) * t226) * t201; 0.2e1 * (-t262 / 0.2e1 + t261 / 0.2e1) * m(6) + 0.2e1 * (t271 * t28 + t272 * t29) * m(5) + 0.2e1 * (t271 * t55 + t272 * t56) * m(4); -t168 * t273 + 0.2e1 * (m(5) * (t160 * t28 + t162 * t29) / 0.2e1 + m(6) * (t10 * t160 + t11 * t162) / 0.2e1) * t166; -t6 * t216 / 0.2e1 + (-t32 * t168 + t178) * t189 + (-t12 * t168 + t173) * t200 + (qJD(5) * t172 + t13 * t150) * t248 / 0.2e1 + (-t31 * t168 + t179) * t188 + (-t13 * t168 + t172) * t202 - t168 * (qJD(5) * t171 + t25) / 0.2e1 + t150 * (-t26 * t168 + t171) / 0.2e1 + ((t106 * t241 + t107 * t242 + t120 * t245) * t150 + (t106 * t274 + t176 * t107 + t162 * t181) * t226) * t201 + ((-t104 * t241 + t105 * t242 + t120 * t248) * t150 + (-t104 * t274 + t105 * t176 + t160 * t181) * t226) * t203 - t150 * (-t168 * t120 * t150 + ((-t159 * t241 + t161 * t242) * t150 + ((-t159 * t274 + t161 * t176) * t166 - t190 * t168) * qJD(5)) * t166) / 0.2e1 + (qJD(2) * t5 + qJD(5) * t173 + t12 * t150) * t245 / 0.2e1 + ((-t10 * t68 + t11 * t66 + t21 * t40 - t22 * t39) * t168 + (t9 * t191 + t27 * (-t231 * t68 - t232 * t66 + t192) + t193 * t114 + (-t262 + t261 + (t160 * t22 + t258) * qJD(2)) * t103) * t166 - (-t21 * t84 + t22 * t85) * t150 - (t27 * (-t160 * t85 + t162 * t84) + t193 * t180) * t226) * m(6);];
tauc = t1(:);
