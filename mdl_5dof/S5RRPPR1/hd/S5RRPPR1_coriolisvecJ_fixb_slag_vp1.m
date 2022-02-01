% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:27
% DurationCPUTime: 3.43s
% Computational Cost: add. (8624->351), mult. (5274->422), div. (0->0), fcn. (3876->10), ass. (0->228)
t171 = qJ(1) + qJ(2);
t164 = pkin(8) + t171;
t159 = cos(t164);
t166 = cos(t171);
t161 = pkin(2) * t166;
t322 = -t159 * rSges(4,1) - t161;
t175 = sin(qJ(1));
t291 = pkin(1) * qJD(1);
t252 = t175 * t291;
t165 = sin(t171);
t122 = rSges(3,1) * t165 + rSges(3,2) * t166;
t170 = qJD(1) + qJD(2);
t279 = t122 * t170;
t93 = -t252 - t279;
t321 = 2 * qJD(5);
t173 = cos(pkin(9));
t160 = pkin(4) * t173 + pkin(3);
t128 = t159 * t160;
t158 = sin(t164);
t174 = -pkin(7) - qJ(4);
t275 = t158 * t174;
t169 = pkin(9) + qJ(5);
t163 = cos(t169);
t273 = t159 * t163;
t212 = rSges(6,1) * t273 + t158 * rSges(6,3);
t162 = sin(t169);
t274 = t159 * t162;
t253 = rSges(6,2) * t274;
t79 = -t253 + t212;
t320 = t128 - t275 + t79;
t155 = t159 * pkin(3);
t315 = t158 * qJ(4) + t155;
t237 = t315 + t161;
t293 = rSges(5,2) * sin(pkin(9));
t254 = t159 * t293;
t295 = rSges(5,1) * t173;
t139 = t159 * t295;
t316 = -t158 * rSges(5,3) - t139;
t319 = t237 - t254 - t316;
t138 = t158 * t293;
t272 = t159 * t170;
t267 = rSges(5,3) * t272 + t170 * t138;
t256 = t158 * t295;
t264 = t159 * rSges(5,3) + t138;
t83 = t256 - t264;
t318 = t170 * t83 + t267;
t317 = -rSges(4,2) * t158 - t322;
t156 = Icges(6,4) * t163;
t215 = -Icges(6,2) * t162 + t156;
t114 = Icges(6,1) * t162 + t156;
t143 = qJD(4) * t158;
t278 = t158 * t162;
t129 = rSges(6,2) * t278;
t292 = rSges(6,2) * t163;
t250 = qJD(5) * t292;
t206 = rSges(6,3) * t272 + t170 * t129 - t159 * t250;
t260 = qJD(5) * t162;
t247 = t159 * t260;
t146 = t159 * qJ(4);
t271 = t159 * t174;
t65 = -t271 - t146 + (pkin(3) - t160) * t158;
t277 = t158 * t163;
t255 = rSges(6,1) * t277;
t266 = t159 * rSges(6,3) + t129;
t78 = t255 - t266;
t67 = t170 * t78;
t314 = -rSges(6,1) * t247 - t170 * t65 + t143 + t206 + t67;
t120 = rSges(6,1) * t162 + t292;
t262 = qJD(5) * t158;
t95 = t120 * t262;
t313 = (t237 - t315 + t320) * t170 - t95;
t111 = Icges(6,5) * t163 - Icges(6,6) * t162;
t110 = Icges(6,5) * t162 + Icges(6,6) * t163;
t192 = Icges(6,3) * t170 - qJD(5) * t110;
t204 = t215 * t159;
t74 = Icges(6,6) * t158 + t204;
t288 = t162 * t74;
t284 = Icges(6,4) * t162;
t115 = Icges(6,1) * t163 - t284;
t205 = t115 * t159;
t76 = Icges(6,5) * t158 + t205;
t219 = -t163 * t76 + t288;
t276 = t158 * t170;
t312 = -t111 * t276 + t159 * t192 + t170 * t219;
t203 = t111 * t159;
t73 = Icges(6,4) * t277 - Icges(6,2) * t278 - Icges(6,6) * t159;
t289 = t162 * t73;
t127 = Icges(6,4) * t278;
t75 = Icges(6,1) * t277 - Icges(6,5) * t159 - t127;
t220 = -t163 * t75 + t289;
t311 = t158 * t192 + (t203 + t220) * t170;
t112 = Icges(6,2) * t163 + t284;
t214 = t112 * t162 - t114 * t163;
t310 = t111 * qJD(5) + t170 * t214;
t71 = Icges(6,5) * t277 - Icges(6,6) * t278 - Icges(6,3) * t159;
t28 = -t158 * t220 - t159 * t71;
t297 = -Icges(6,2) * t277 - t127 + t75;
t299 = t114 * t158 + t73;
t309 = -t162 * t297 - t163 * t299;
t168 = t170 ^ 2;
t308 = t158 / 0.2e1;
t307 = -t159 / 0.2e1;
t306 = t170 / 0.2e1;
t305 = pkin(1) * t175;
t304 = pkin(2) * t165;
t303 = pkin(2) * t168;
t302 = pkin(3) * t158;
t176 = cos(qJ(1));
t167 = t176 * pkin(1);
t301 = -t158 * t71 - t75 * t273;
t72 = Icges(6,3) * t158 + t203;
t300 = t158 * t72 + t76 * t273;
t298 = -t114 * t159 - t74;
t296 = -t112 * t159 + t76;
t294 = rSges(6,1) * t163;
t157 = t166 * rSges(3,1);
t281 = t110 * t159;
t47 = -t158 * t214 - t281;
t287 = t47 * t170;
t104 = -t146 + t302;
t96 = t170 * t104;
t286 = t143 - t96;
t282 = t110 * t158;
t280 = t111 * t170;
t270 = t165 * t170;
t269 = -t112 + t115;
t268 = t114 + t215;
t134 = qJ(4) * t272;
t265 = t134 + t143;
t263 = qJD(4) * t170;
t261 = qJD(5) * t159;
t144 = qJD(4) * t159;
t251 = t176 * t291;
t228 = -t144 + t251;
t45 = t170 * t319 + t228;
t259 = t45 * t304;
t177 = qJD(1) ^ 2;
t258 = t177 * t305;
t257 = t177 * t167;
t249 = -t170 * t253 + (-rSges(6,1) * t260 - t250) * t158;
t248 = t120 * t261;
t244 = -pkin(3) - t295;
t243 = -t262 / 0.2e1;
t240 = t261 / 0.2e1;
t238 = -t104 - t304;
t105 = rSges(4,1) * t158 + rSges(4,2) * t159;
t201 = -t105 - t304;
t57 = t76 * t277;
t235 = t159 * t72 - t57;
t234 = -t71 + t288;
t123 = -rSges(3,2) * t165 + t157;
t99 = -rSges(3,2) * t270 + t157 * t170;
t229 = t143 - t252;
t226 = t143 - t248;
t224 = -rSges(6,2) * t162 + t294;
t26 = -t252 + (t238 + t65 - t78) * t170 + t226;
t27 = t228 + t313;
t222 = -t158 * t27 - t159 * t26;
t221 = t158 * t78 + t159 * t79;
t42 = t162 * t75 + t163 * t73;
t43 = t162 * t76 + t163 * t74;
t124 = t170 * t275;
t218 = t124 + t144 - t249;
t213 = -t158 * t160 - t271;
t29 = -t278 * t74 - t235;
t210 = (t158 * t29 - t159 * t28) * qJD(5);
t30 = -t274 * t73 - t301;
t31 = -t274 * t74 + t300;
t209 = (t158 * t31 - t159 * t30) * qJD(5);
t208 = -t165 * t303 - t258;
t207 = -t166 * t303 - t257;
t202 = -pkin(2) * t270 - t252;
t200 = t159 * t263 + t207;
t199 = -t162 * t296 + t163 * t298;
t198 = t158 * t263 + t208 + t170 * (-pkin(3) * t276 + t265);
t196 = t202 + t286;
t195 = (-t162 * t268 + t163 * t269) * t170;
t194 = Icges(6,5) * t170 - qJD(5) * t114;
t193 = Icges(6,6) * t170 - qJD(5) * t112;
t191 = t161 + t320;
t190 = t158 * t244 + t146 + t264 - t304;
t55 = (-t163 * t276 - t247) * rSges(6,1) + t206;
t56 = t170 * t212 + t249;
t189 = (t55 + t67) * t159 + (-t170 * t79 + t56) * t158;
t48 = -t159 * t214 + t282;
t46 = t48 * t170;
t10 = t46 + t209;
t101 = t215 * qJD(5);
t102 = t115 * qJD(5);
t52 = t158 * t193 + t170 * t204;
t54 = t158 * t194 + t170 * t205;
t16 = -qJD(5) * t220 + t162 * t54 + t163 * t52;
t51 = t159 * t193 - t215 * t276;
t53 = -t115 * t276 + t159 * t194;
t17 = -qJD(5) * t219 + t162 * t53 + t163 * t51;
t181 = -t101 * t162 + t102 * t163 + t110 * t170 + (-t112 * t163 - t114 * t162) * qJD(5);
t20 = t310 * t158 + t181 * t159;
t21 = t181 * t158 - t310 * t159;
t9 = t210 + t287;
t188 = (t46 + ((t29 - t57 + (t72 + t289) * t159 + t301) * t159 + t300 * t158) * qJD(5)) * t240 + (-qJD(5) * t214 + t101 * t163 + t102 * t162) * t170 + (-t287 + ((t159 * t234 - t300 + t31) * t159 + (t158 * t234 + t235 + t30) * t158) * qJD(5) + t9) * t243 + (t17 + t20) * t262 / 0.2e1 - (t16 + t21 + t10) * t261 / 0.2e1 + ((t42 + t47) * t158 + (t43 + t48) * t159) * qJD(5) * t306;
t184 = -t304 - t271 + (-t160 - t294) * t158 + t266;
t183 = -qJD(5) * t43 - t162 * t51 + t163 * t53 + t170 * t72;
t182 = -qJD(5) * t42 - t162 * t52 + t163 * t54 + t170 * t71;
t68 = t170 * t201 - t252;
t69 = t170 * t317 + t251;
t180 = (t69 * t201 + t68 * t322) * t170;
t44 = (t238 - t83) * t170 + t229;
t179 = (t44 * (-t139 - t155 - t161) - t259 + (t44 * (-rSges(5,3) - qJ(4)) + t45 * t244) * t158) * t170;
t178 = (t26 * (-t128 - t212 - t161) + t27 * (t213 - t255 - t304)) * t170;
t135 = rSges(4,2) * t276;
t117 = t170 * t254;
t103 = t224 * qJD(5);
t97 = t170 * t105;
t94 = t123 * t170 + t251;
t92 = t120 * t159;
t91 = t120 * t158;
t82 = -t170 * t99 - t257;
t81 = -t170 * t279 - t258;
t80 = t170 * t315 - t144;
t61 = -t170 * (rSges(4,1) * t272 - t135) + t207;
t60 = -t105 * t168 + t208;
t41 = qJD(5) * t221 + qJD(3);
t33 = (t316 * t170 + t117 - t80) * t170 + t200;
t32 = t170 * (-t170 * t256 + t267) + t198;
t15 = -t103 * t261 + (t95 + t124 - t56 - t80 + (t315 - t128) * t170) * t170 + t200;
t14 = -t103 * t262 + (-t248 - t134 + t55 + (t213 + t302) * t170) * t170 + t198;
t11 = t189 * qJD(5);
t1 = [m(3) * (t82 * (-t122 - t305) + t81 * (t123 + t167) + (-t99 - t251 + t94) * t93) + t188 + (t15 * (t184 - t305) + t26 * (t218 - t251) + t14 * (t167 + t191) + t178 + (-t196 + t26 + t248 - t252 + t314) * t27) * m(6) + (t33 * (t190 - t305) + t44 * (t117 - t228) + t32 * (t167 + t319) + t179 + (-t196 + t44 + t134 + t229 + t318) * t45) * m(5) + (t61 * (t201 - t305) + t68 * (t135 - t251) + t60 * (t167 + t317) + t180 + (-t202 + t68 + t97 - t252) * t69) * m(4); t188 + (t14 * t191 + t15 * t184 + t178 + (t304 * t170 - t226 + t314 + t96) * t27 + (-t144 + t218 + t313) * t26) * m(6) + (-(-t319 * t44 - t259) * t170 + t33 * t190 + t44 * t117 + t32 * t319 + t179 + (-t286 + t265 + t318) * t45) * m(5) + (t69 * t97 - (-t304 * t69 - t317 * t68) * t170 + t68 * t135 + t61 * t201 + t60 * t317 + t180) * m(4) + (-(-t122 * t94 - t123 * t93) * t170 - t122 * t82 + t123 * t81 - t93 * t99 - t94 * t279) * m(3); m(6) * t11; 0.2e1 * (t14 * t307 + t15 * t308) * m(6) + 0.2e1 * (t307 * t32 + t308 * t33) * m(5); ((t170 * t43 - t16) * t159 + (t170 * t42 + t17) * t158) * t306 + ((-t262 * t281 + t280) * t158 + (t195 + (-t309 * t159 + (t282 + t199) * t158) * qJD(5)) * t159) * t243 + ((-t261 * t282 - t280) * t159 + (t195 + (t199 * t158 + (-t309 + t281) * t159) * qJD(5)) * t158) * t240 - t170 * ((t162 * t269 + t163 * t268) * t170 + ((t158 * t296 - t159 * t297) * t163 + (t158 * t298 + t159 * t299) * t162) * qJD(5)) / 0.2e1 + (t170 * t20 + ((-t311 * t158 - t182 * t159 + t170 * t31) * t159 + (t312 * t158 + t183 * t159 + t170 * t30) * t158) * t321) * t308 + (t170 * t21 + ((-t182 * t158 + t311 * t159 + t170 * t29) * t159 + (t183 * t158 - t312 * t159 + t170 * t28) * t158) * t321) * t307 + (t9 + t210) * t276 / 0.2e1 + (t10 + t209) * t272 / 0.2e1 + (t11 * t221 + t41 * t189 + t222 * t103 + ((-t170 * t27 - t15) * t159 + (t170 * t26 - t14) * t158) * t120 - (t26 * t91 - t27 * t92) * t170 - (t41 * (-t158 * t91 - t159 * t92) + t222 * t224) * qJD(5)) * m(6);];
tauc = t1(:);
