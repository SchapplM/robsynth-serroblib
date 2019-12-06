% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:17
% DurationCPUTime: 8.79s
% Computational Cost: add. (8798->444), mult. (23999->596), div. (0->0), fcn. (27434->10), ass. (0->217)
t198 = sin(pkin(9));
t201 = cos(pkin(9));
t202 = cos(pkin(7));
t200 = sin(pkin(7));
t297 = cos(pkin(8));
t254 = t200 * t297;
t167 = t198 * t254 + t202 * t201;
t142 = qJD(5) * t167 + qJD(1);
t204 = sin(qJ(1));
t253 = t204 * t297;
t199 = sin(pkin(8));
t206 = cos(qJ(1));
t287 = t206 * t199;
t170 = -t202 * t253 + t287;
t290 = t200 * t204;
t134 = t170 * t198 + t201 * t290;
t252 = t206 * t297;
t288 = t204 * t199;
t172 = t202 * t252 + t288;
t289 = t200 * t206;
t136 = t172 * t198 - t201 * t289;
t135 = t170 * t201 - t198 * t290;
t169 = t202 * t288 + t252;
t203 = sin(qJ(5));
t205 = cos(qJ(5));
t102 = -t135 * t203 - t169 * t205;
t103 = t135 * t205 - t169 * t203;
t51 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t134;
t296 = Icges(6,4) * t103;
t54 = Icges(6,2) * t102 + Icges(6,6) * t134 + t296;
t98 = Icges(6,4) * t102;
t57 = Icges(6,1) * t103 + Icges(6,5) * t134 + t98;
t14 = t102 * t54 + t103 * t57 + t134 * t51;
t137 = t172 * t201 + t198 * t289;
t171 = t202 * t287 - t253;
t105 = t137 * t205 + t171 * t203;
t106 = t137 * t203 - t171 * t205;
t52 = Icges(6,5) * t105 - Icges(6,6) * t106 + Icges(6,3) * t136;
t295 = Icges(6,4) * t105;
t56 = Icges(6,2) * t106 - Icges(6,6) * t136 - t295;
t99 = Icges(6,4) * t106;
t58 = Icges(6,1) * t105 + Icges(6,5) * t136 - t99;
t15 = -t102 * t56 + t103 * t58 + t134 * t52;
t235 = t134 * t14 + t136 * t15;
t168 = -t202 * t198 + t201 * t254;
t291 = t199 * t200;
t132 = -t168 * t203 + t205 * t291;
t133 = t168 * t205 + t203 * t291;
t82 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t167;
t294 = Icges(6,4) * t133;
t83 = Icges(6,2) * t132 + Icges(6,6) * t167 + t294;
t127 = Icges(6,4) * t132;
t84 = Icges(6,1) * t133 + Icges(6,5) * t167 + t127;
t26 = t102 * t83 + t103 * t84 + t134 * t82;
t5 = qJD(5) * t235 + t26 * t142;
t272 = qJD(5) * t136;
t256 = t272 / 0.2e1;
t16 = t105 * t57 - t106 * t54 + t136 * t51;
t17 = t105 * t58 + t106 * t56 + t52 * t136;
t234 = t134 * t16 + t136 * t17;
t21 = -t132 * t56 + t133 * t58 + t167 * t52;
t85 = rSges(6,1) * t133 + rSges(6,2) * t132 + rSges(6,3) * t167;
t328 = t136 * t85;
t61 = rSges(6,1) * t105 - rSges(6,2) * t106 + rSges(6,3) * t136;
t327 = t142 * t61;
t27 = t105 * t84 - t106 * t83 + t136 * t82;
t320 = t134 * t85;
t275 = qJD(3) * t206;
t186 = t200 * t275;
t307 = pkin(2) * t202;
t238 = qJ(3) * t200 + t307;
t224 = qJD(1) * t238;
t193 = qJD(2) * t204;
t195 = t206 * qJ(2);
t282 = qJD(1) * (-pkin(1) * t204 + t195) + t193;
t250 = -t204 * t224 + t186 + t282;
t251 = t170 * pkin(3) - qJ(4) * t169;
t274 = qJD(4) * t171;
t311 = qJD(1) * t251 + t250 + t274;
t310 = qJ(2) * qJD(1);
t248 = -pkin(4) * t137 - pkin(6) * t136;
t221 = -rSges(5,1) * t137 + rSges(5,2) * t136 - rSges(5,3) * t171;
t116 = pkin(3) * t172 + t171 * qJ(4);
t60 = t103 * rSges(6,1) + t102 * rSges(6,2) + t134 * rSges(6,3);
t95 = t135 * pkin(4) + pkin(6) * t134;
t309 = -qJD(1) * t95 - t142 * t60 - t311;
t86 = t135 * rSges(5,1) - t134 * rSges(5,2) - t169 * rSges(5,3);
t46 = qJD(1) * t86 + t311;
t308 = pkin(1) * t206;
t304 = qJD(5) / 0.2e1;
t303 = Icges(6,1) * t132 - t294 - t83;
t302 = -Icges(6,2) * t133 + t127 + t84;
t301 = rSges(3,1) * t202;
t299 = -rSges(3,3) - qJ(2);
t298 = -rSges(4,3) - qJ(3);
t153 = t171 * qJD(1);
t154 = t172 * qJD(1);
t286 = -t154 * rSges(4,1) + t153 * rSges(4,2);
t194 = qJD(2) * t206;
t180 = qJ(2) * t204 + t308;
t177 = qJD(1) * t180;
t279 = t194 - t177;
t285 = (t194 + t279) * qJD(1);
t284 = t170 * rSges(4,1) + t169 * rSges(4,2);
t174 = t238 * t206;
t283 = -t174 - t180;
t281 = rSges(3,2) * t290 + t206 * rSges(3,3);
t278 = qJD(1) * t204;
t280 = -pkin(1) * t278 + t193;
t277 = qJD(1) * t206;
t276 = qJD(3) * t204;
t273 = qJD(5) * t134;
t157 = t169 * qJD(4);
t264 = t200 * t277;
t123 = -t154 * t198 + t201 * t264;
t124 = -t154 * t201 - t198 * t264;
t69 = -qJD(5) * t103 - t124 * t203 - t153 * t205;
t70 = qJD(5) * t102 + t124 * t205 - t153 * t203;
t36 = t70 * rSges(6,1) + t69 * rSges(6,2) + t123 * rSges(6,3);
t271 = t206 * t301;
t269 = -t116 + t283;
t268 = t124 * rSges(5,1) - t123 * rSges(5,2) - t153 * rSges(5,3);
t263 = t200 * t276;
t267 = qJD(1) * (-t206 * t224 - t263) + t285;
t265 = t200 * t278;
t141 = -qJ(3) * t265 - t278 * t307 + t186;
t266 = -pkin(1) - t307;
t262 = -pkin(1) - t301;
t152 = t170 * qJD(1);
t121 = t152 * t198 + t201 * t265;
t261 = t121 * t304;
t260 = t123 * t304;
t258 = t273 / 0.2e1;
t257 = -t272 / 0.2e1;
t255 = t124 * pkin(4) + pkin(6) * t123;
t151 = t169 * qJD(1);
t97 = pkin(3) * t152 - t151 * qJ(4) + t274;
t122 = t152 * t201 - t198 * t265;
t249 = -pkin(4) * t122 - pkin(6) * t121;
t247 = t195 + t251;
t246 = t194 - t263;
t166 = qJ(2) * t277 + t280;
t245 = -t141 - t166 - t193;
t225 = -t154 * pkin(3) - qJ(4) * t153 - t157;
t244 = qJD(1) * t225 - qJD(4) * t151 + t267;
t241 = -rSges(3,2) * t200 + t301;
t240 = -rSges(4,1) * t152 - rSges(4,2) * t151;
t239 = -rSges(4,1) * t172 + rSges(4,2) * t171;
t125 = t132 * qJD(5);
t126 = t133 * qJD(5);
t90 = rSges(6,1) * t125 - rSges(6,2) * t126;
t237 = t121 * t85 + t136 * t90;
t236 = -t123 * t85 - t134 * t90;
t233 = t134 * t61 - t136 * t60;
t232 = t134 * (Icges(6,5) * t102 - Icges(6,6) * t103) + t136 * (-Icges(6,5) * t106 - Icges(6,6) * t105);
t230 = -t141 - t280;
t229 = -qJD(3) * t200 - t310;
t228 = -pkin(1) - t238;
t227 = -t157 + t246;
t226 = -rSges(3,3) * t204 - t271;
t223 = -qJD(1) * t174 - t177 + t246;
t67 = -qJD(5) * t105 - t122 * t203 - t151 * t205;
t68 = -qJD(5) * t106 + t122 * t205 - t151 * t203;
t35 = rSges(6,1) * t68 + rSges(6,2) * t67 + rSges(6,3) * t121;
t222 = -rSges(5,1) * t122 + rSges(5,2) * t121 + rSges(5,3) * t151;
t220 = t194 + t225;
t218 = (Icges(6,1) * t102 - t296 - t54) * t134 + (-Icges(6,1) * t106 - t295 + t56) * t136;
t217 = (-Icges(6,2) * t103 + t57 + t98) * t134 + (-Icges(6,2) * t105 + t58 - t99) * t136;
t216 = -rSges(4,3) * t289 + t239;
t215 = -qJD(1) * t116 - t157 + t223;
t214 = t245 - t97 - t186;
t213 = t230 - t97;
t29 = Icges(6,5) * t68 + Icges(6,6) * t67 + Icges(6,3) * t121;
t30 = Icges(6,5) * t70 + Icges(6,6) * t69 + Icges(6,3) * t123;
t31 = Icges(6,4) * t68 + Icges(6,2) * t67 + Icges(6,6) * t121;
t32 = Icges(6,4) * t70 + Icges(6,2) * t69 + Icges(6,6) * t123;
t33 = Icges(6,1) * t68 + Icges(6,4) * t67 + Icges(6,5) * t121;
t34 = Icges(6,1) * t70 + Icges(6,4) * t69 + Icges(6,5) * t123;
t212 = (t105 * t34 - t106 * t32 + t121 * t51 + t136 * t30 + t54 * t67 + t57 * t68) * t134 + t121 * t17 + t123 * t16 + t136 * (t105 * t33 - t106 * t31 + t121 * t52 + t136 * t29 - t56 * t67 + t58 * t68);
t211 = t121 * t15 + t123 * t14 + t134 * (t102 * t32 + t103 * t34 + t123 * t51 + t134 * t30 + t54 * t69 + t57 * t70) + t136 * (t102 * t31 + t103 * t33 + t123 * t52 + t134 * t29 - t56 * t69 + t58 * t70);
t20 = t132 * t54 + t133 * t57 + t167 * t51;
t7 = t125 * t57 - t126 * t54 + t132 * t32 + t133 * t34 + t167 * t30;
t8 = t125 * t58 + t126 * t56 + t132 * t31 + t133 * t33 + t167 * t29;
t210 = t121 * t21 + t123 * t20 + t134 * t7 + t136 * t8;
t209 = -t121 * t60 + t123 * t61 + t134 * t35 - t136 * t36;
t189 = rSges(3,2) * t289;
t184 = rSges(3,2) * t264;
t150 = -t189 - t226;
t144 = qJD(4) * t153;
t112 = t194 + (-t150 - t180) * qJD(1);
t111 = qJD(1) * (-t204 * t301 + t281) + t282;
t109 = (qJD(1) * t226 + t184) * qJD(1) + t285;
t108 = (-rSges(3,3) * t277 - t166 + (qJD(1) * t241 - qJD(2)) * t204) * qJD(1);
t91 = Icges(6,5) * t132 - Icges(6,6) * t133;
t89 = Icges(6,1) * t125 - Icges(6,4) * t126;
t88 = Icges(6,4) * t125 - Icges(6,2) * t126;
t87 = Icges(6,5) * t125 - Icges(6,6) * t126;
t81 = (t216 + t283) * qJD(1) + t246;
t80 = qJD(1) * (-rSges(4,3) * t290 + t284) + t250;
t78 = -rSges(6,1) * t106 - rSges(6,2) * t105;
t77 = rSges(6,1) * t102 - rSges(6,2) * t103;
t66 = ((-rSges(4,3) * t277 - t276) * t200 + t286) * qJD(1) + t267;
t65 = ((rSges(4,3) * t278 - t275) * t200 + t240 + t245) * qJD(1);
t47 = (t221 + t269) * qJD(1) + t227;
t38 = (-t263 + t268) * qJD(1) + t244;
t37 = -t144 + (t214 + t222) * qJD(1);
t28 = -qJD(3) * t202 + qJD(4) * t291 + qJD(5) * t233;
t24 = t85 * t272 - t327 + (t248 + t269) * qJD(1) + t227;
t23 = -t273 * t85 - t309;
t19 = t125 * t84 - t126 * t83 + t132 * t88 + t133 * t89 + t167 * t87;
t18 = t19 * t142;
t13 = -t142 * t35 - t144 + t237 * qJD(5) + (t214 + t249) * qJD(1);
t12 = t142 * t36 + t236 * qJD(5) + (t255 - t263) * qJD(1) + t244;
t11 = t102 * t88 + t103 * t89 + t123 * t82 + t134 * t87 + t69 * t83 + t70 * t84;
t10 = t105 * t89 - t106 * t88 + t121 * t82 + t136 * t87 + t67 * t83 + t68 * t84;
t9 = t209 * qJD(5);
t1 = [t5 * t257 + t18 + (t21 + t27) * t261 + (t20 + t26) * t260 + (t7 + t11) * t258 + (t13 * (-t116 + t248 - t61) + t24 * (t213 + t249 - t35) + t12 * (t247 + t95 + t60) + t23 * (t220 + t255 + t36) + (-t13 * qJ(2) + t12 * t228 + t229 * t23) * t204 + (t13 * t228 + (-t24 * qJ(2) + t228 * t23) * qJD(1)) * t206 - t24 * t309 - t23 * (qJD(1) * t248 + t215 - t327) - (t23 * t328 + t24 * t320) * qJD(5)) * m(6) + (t37 * (-t116 + t221) + t47 * (t213 + t222) + t38 * (t247 + t86) + (-t37 * qJ(2) + t228 * t38) * t204 + (t37 * t228 - t47 * t310) * t206 + (t229 * t204 - t215 + t220 + t268 + t47 + (t228 * t206 - t221) * qJD(1)) * t46) * m(5) + (t65 * t239 + t81 * (t230 + t240) + t66 * (t195 + t284) + t80 * (t194 + t286) + (t66 * t266 + (-qJD(1) * t80 - t65) * qJ(2) + (t81 * rSges(4,3) * qJD(1) - t80 * qJD(3) + t298 * t66) * t200) * t204 + (t65 * (t200 * t298 + t266) + (-t81 * qJ(2) + t80 * (-rSges(4,3) * t200 + t228)) * qJD(1)) * t206 - (qJD(1) * t216 + t223 - t81) * t80) * m(4) + (t108 * (t189 - t271 - t308) - t112 * t280 + t109 * (t195 + t281) + t111 * (t184 + t194) + (t108 * t299 + t109 * t262) * t204 + ((t111 * t262 + t112 * t299) * t206 + (t111 * t299 + t112 * t241) * t204) * qJD(1) - (-qJD(1) * t150 - t112 + t279) * t111) * m(3) + (t5 + t8 + t10) * t256; m(3) * (t108 * t206 + t109 * t204) + m(4) * (t204 * t66 + t206 * t65) + m(5) * (t204 * t38 + t206 * t37) + m(6) * (t12 * t204 + t13 * t206); -m(6) * t202 * t9 + 0.2e1 * (m(4) * (-t204 * t65 + t206 * t66) / 0.2e1 + m(5) * (-t204 * t37 + t206 * t38) / 0.2e1 + m(6) * (t12 * t206 - t13 * t204) / 0.2e1) * t200; m(5) * (-t151 * t46 - t153 * t47 - t169 * t37 + t171 * t38) + m(6) * (t12 * t171 - t13 * t169 - t151 * t23 - t153 * t24 + t9 * t291) + 0.2e1 * (-m(5) * (-t169 * t46 - t171 * t47) / 0.2e1 - m(6) * (-t169 * t23 - t171 * t24) / 0.2e1) * qJD(1); t167 * (t210 * qJD(5) + t18) / 0.2e1 + t142 * (t167 * t19 + t210) / 0.2e1 + t123 * t5 / 0.2e1 + t134 * (t211 * qJD(5) + t11 * t142) / 0.2e1 + (t167 * t26 + t235) * t260 + (t11 * t167 + t211) * t258 + t121 * (qJD(5) * t234 + t142 * t27) / 0.2e1 + t136 * (t212 * qJD(5) + t10 * t142) / 0.2e1 + (t167 * t27 + t234) * t261 + (t10 * t167 + t212) * t256 - t142 * ((t302 * t132 + t303 * t133 + t167 * t91) * t142 + (t132 * t217 + t133 * t218 + t167 * t232) * qJD(5)) / 0.2e1 - ((t102 * t302 + t103 * t303 + t134 * t91) * t142 + (t102 * t217 + t103 * t218 + t134 * t232) * qJD(5)) * t273 / 0.2e1 + ((t105 * t303 - t106 * t302 + t136 * t91) * t142 + (t105 * t218 - t106 * t217 + t136 * t232) * qJD(5)) * t257 + (t9 * t233 + t28 * t209 + t13 * (-t167 * t61 + t328) + t24 * (-t167 * t35 + t237) + t12 * (t167 * t60 - t320) + t23 * (t167 * t36 + t236) - (t23 * t77 - t24 * t78) * t142 - (t28 * (t134 * t78 - t136 * t77) + (-t134 * t23 + t136 * t24) * (rSges(6,1) * t132 - rSges(6,2) * t133)) * qJD(5)) * m(6);];
tauc = t1(:);
