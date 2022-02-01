% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:38
% DurationCPUTime: 5.20s
% Computational Cost: add. (11071->389), mult. (9267->525), div. (0->0), fcn. (8721->10), ass. (0->216)
t180 = sin(pkin(9));
t181 = cos(pkin(9));
t182 = sin(qJ(5));
t184 = cos(qJ(5));
t179 = qJ(1) + qJ(2);
t173 = pkin(8) + t179;
t170 = sin(t173);
t171 = cos(t173);
t266 = t181 * t182;
t115 = t170 * t266 + t171 * t184;
t265 = t181 * t184;
t116 = t170 * t265 - t171 * t182;
t272 = t170 * t180;
t68 = Icges(6,5) * t116 - Icges(6,6) * t115 + Icges(6,3) * t272;
t109 = Icges(6,4) * t116;
t71 = -Icges(6,2) * t115 + Icges(6,6) * t272 + t109;
t108 = Icges(6,4) * t115;
t75 = -Icges(6,1) * t116 - Icges(6,5) * t272 + t108;
t30 = -(t182 * t71 + t184 * t75) * t180 - t181 * t68;
t212 = (-rSges(6,1) * t182 - rSges(6,2) * t184) * t180;
t132 = qJD(5) * t212;
t257 = qJD(5) * t180;
t312 = t132 * t257;
t311 = -t115 * t71 - t116 * t75;
t117 = t170 * t184 - t171 * t266;
t118 = t170 * t182 + t171 * t265;
t310 = t117 * t71 - t118 * t75;
t270 = t171 * t180;
t307 = t68 * t270;
t168 = t171 * pkin(3);
t135 = t170 * qJ(4) + t168;
t175 = cos(t179);
t172 = pkin(2) * t175;
t297 = t172 + t135;
t269 = t171 * t181;
t153 = rSges(5,1) * t269;
t298 = -t170 * rSges(5,3) - t153;
t306 = -rSges(5,2) * t270 + t297 - t298;
t305 = pkin(4) * t269 + pkin(7) * t270 + t297;
t304 = -t171 * rSges(4,1) - t172;
t183 = sin(qJ(1));
t282 = pkin(1) * qJD(1);
t252 = t183 * t282;
t174 = sin(t179);
t143 = rSges(3,1) * t174 + rSges(3,2) * t175;
t178 = qJD(1) + qJD(2);
t274 = t143 * t178;
t113 = -t252 - t274;
t160 = -qJD(5) * t181 + t178;
t128 = -rSges(6,3) * t181 + (rSges(6,1) * t184 - rSges(6,2) * t182) * t180;
t245 = t128 * t257;
t223 = rSges(6,1) * t116 - rSges(6,2) * t115;
t77 = rSges(6,3) * t272 + t223;
t303 = -t160 * t77 + t170 * t245;
t287 = pkin(4) * t181;
t206 = -t287 - pkin(3) + (-rSges(6,3) - pkin(7)) * t180;
t225 = pkin(7) * t180 + t287;
t162 = qJD(4) * t171;
t100 = t135 * t178 - t162;
t195 = -t225 * t178 + t245;
t185 = cos(qJ(1));
t176 = t185 * pkin(1);
t186 = qJD(1) ^ 2;
t254 = t186 * t176;
t177 = t178 ^ 2;
t288 = pkin(2) * t177;
t207 = -t175 * t288 - t254;
t258 = qJD(4) * t178;
t203 = t171 * t258 + t207;
t84 = -t116 * qJD(5) + t117 * t178;
t85 = -t115 * qJD(5) + t118 * t178;
t226 = rSges(6,1) * t85 + rSges(6,2) * t84;
t267 = t178 * t180;
t248 = t171 * t267;
t49 = rSges(6,3) * t248 + t226;
t22 = t170 * t312 - t160 * t49 + (t195 * t171 - t100) * t178 + t203;
t278 = t22 * t170;
t251 = t185 * t282;
t79 = t118 * rSges(6,1) + t117 * rSges(6,2) + rSges(6,3) * t270;
t294 = t160 * t79 - t171 * t245 + t305 * t178 - t162;
t35 = t251 + t294;
t279 = t174 * t35;
t119 = t225 * t170;
t161 = qJD(4) * t170;
t230 = t161 - t252;
t164 = t171 * qJ(4);
t133 = pkin(3) * t170 - t164;
t289 = pkin(2) * t174;
t240 = -t133 - t289;
t34 = (-t119 + t240) * t178 + t230 + t303;
t302 = ((-t175 * t34 - t279) * pkin(2) + t34 * t206 * t171 + (-t34 * qJ(4) + t35 * (-rSges(6,3) * t180 - pkin(3) - t225)) * t170) * t178 + t206 * t278;
t26 = t307 + t310;
t301 = t26 - t307;
t283 = rSges(5,1) * t181;
t253 = t170 * t283;
t260 = rSges(5,2) * t272 + t171 * rSges(5,3);
t105 = t253 - t260;
t249 = t170 * t267;
t271 = t171 * t178;
t262 = rSges(5,2) * t249 + rSges(5,3) * t271;
t300 = t178 * t105 + t262;
t299 = -rSges(4,2) * t170 - t304;
t148 = qJ(4) * t271;
t268 = t174 * t178;
t205 = -pkin(2) * t268 - t252;
t259 = -t178 * t133 + t161;
t296 = t148 + t230 - t205 - t259;
t261 = t148 + t161;
t295 = t261 - t259;
t70 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t270;
t277 = Icges(6,4) * t118;
t73 = Icges(6,2) * t117 + Icges(6,6) * t270 + t277;
t110 = Icges(6,4) * t117;
t76 = Icges(6,1) * t118 + Icges(6,5) * t270 + t110;
t25 = -t115 * t73 + t116 * t76 + t70 * t272;
t82 = -t118 * qJD(5) + t115 * t178;
t83 = t117 * qJD(5) - t116 * t178;
t284 = t83 * rSges(6,1) + t82 * rSges(6,2);
t293 = t178 * t119 + t284 - t303;
t292 = (-Icges(6,2) * t116 - t108 - t75) * t170 + (-Icges(6,2) * t118 + t110 + t76) * t171;
t291 = -t171 / 0.2e1;
t290 = pkin(1) * t183;
t121 = -Icges(6,3) * t181 + (Icges(6,5) * t184 - Icges(6,6) * t182) * t180;
t275 = Icges(6,4) * t184;
t122 = -Icges(6,6) * t181 + (-Icges(6,2) * t182 + t275) * t180;
t276 = Icges(6,4) * t182;
t123 = -Icges(6,5) * t181 + (Icges(6,1) * t184 - t276) * t180;
t50 = -t115 * t122 + t116 * t123 + t121 * t272;
t281 = t160 * t50;
t24 = t68 * t272 + t311;
t280 = t170 * t24;
t169 = t175 * rSges(3,1);
t273 = t170 * t178;
t139 = (-Icges(6,1) * t182 - t275) * t180;
t264 = -t122 + t139;
t138 = (-Icges(6,2) * t184 - t276) * t180;
t263 = t123 + t138;
t229 = -t162 + t251;
t63 = t178 * t306 + t229;
t256 = t63 * t289;
t27 = t117 * t73 + t118 * t76 + t70 * t270;
t255 = t186 * t290;
t243 = -pkin(3) - t283;
t242 = -t257 / 0.2e1;
t241 = t257 / 0.2e1;
t134 = rSges(4,1) * t170 + rSges(4,2) * t171;
t204 = -t134 - t289;
t237 = t164 - t289;
t144 = -rSges(3,2) * t174 + t169;
t236 = t170 * t242;
t235 = t170 * t241;
t234 = t171 * t242;
t233 = t171 * t241;
t127 = -rSges(3,2) * t268 + t178 * t169;
t222 = t170 * t34 - t171 * t35;
t221 = -t170 * t79 + t171 * t77;
t220 = t170 * (-Icges(6,5) * t115 - Icges(6,6) * t116) + t171 * (Icges(6,5) * t117 - Icges(6,6) * t118);
t219 = t178 * t236;
t218 = t178 * t233;
t216 = t162 - t226;
t213 = t180 * t220;
t211 = (t171 * t25 + t280) * t180;
t210 = (t170 * t26 + t171 * t27) * t180;
t208 = -t174 * t288 - t255;
t137 = (-Icges(6,5) * t182 - Icges(6,6) * t184) * t180;
t202 = (Icges(6,1) * t117 - t277 - t73) * t171 + (-Icges(6,1) * t115 - t109 - t71) * t170;
t200 = t79 + t305;
t199 = -t223 + t237;
t198 = t170 * t258 + t208 + t178 * (-pkin(3) * t273 + t261);
t194 = t243 * t170 + t237 + t260;
t51 = t117 * t122 + t118 * t123 + t121 * t270;
t41 = t51 * t160;
t10 = qJD(5) * t210 + t41;
t43 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t248;
t45 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t248;
t47 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t248;
t13 = -t181 * t43 + (-t182 * t45 + t184 * t47 + (t182 * t75 - t184 * t71) * qJD(5)) * t180;
t42 = Icges(6,5) * t83 + Icges(6,6) * t82 - Icges(6,3) * t249;
t44 = Icges(6,4) * t83 + Icges(6,2) * t82 - Icges(6,6) * t249;
t46 = Icges(6,1) * t83 + Icges(6,4) * t82 - Icges(6,5) * t249;
t14 = -t181 * t42 + (-t182 * t44 + t184 * t46 + (-t182 * t76 - t184 * t73) * qJD(5)) * t180;
t129 = qJD(5) * t137;
t130 = qJD(5) * t138;
t131 = qJD(5) * t139;
t19 = t117 * t130 + t118 * t131 + t122 * t82 + t123 * t83 + (-t121 * t273 + t129 * t171) * t180;
t20 = -t115 * t130 + t116 * t131 + t122 * t84 + t123 * t85 + (t121 * t271 + t129 * t170) * t180;
t31 = -t181 * t70 + (-t182 * t73 + t184 * t76) * t180;
t40 = -t129 * t181 + (-t130 * t182 + t131 * t184 + (-t122 * t184 - t123 * t182) * qJD(5)) * t180;
t37 = t40 * t160;
t9 = qJD(5) * t211 + t281;
t193 = (t41 + ((t24 + t27 - t311) * t171 + t301 * t170) * t257) * t236 + t37 + (t14 + t19) * t233 + (t31 + t51) * t219 + (t30 + t50) * t218 + (t10 + t20 + t13) * t235 + (t9 - t281 + ((-t25 + t301 - t310) * t171 - t280) * t257) * t234;
t192 = ((t178 * t24 - t115 * t44 + t116 * t46 + t73 * t84 + t76 * t85 + (t170 * t42 + t70 * t271) * t180) * t171 + (-t178 * t25 - t115 * t45 + t116 * t47 + t71 * t84 - t75 * t85 + (t170 * t43 + t68 * t271) * t180) * t170) * t180;
t191 = ((t178 * t26 + t117 * t44 + t118 * t46 + t73 * t82 + t76 * t83 + (t171 * t42 - t70 * t273) * t180) * t171 + (-t178 * t27 + t117 * t45 + t118 * t47 + t71 * t82 - t75 * t83 + (t171 * t43 - t68 * t273) * t180) * t170) * t180;
t190 = ((t178 * t30 + t14) * t171 + (-t178 * t31 + t13) * t170) * t180;
t96 = t204 * t178 - t252;
t97 = t178 * t299 + t251;
t189 = (t97 * t204 + t96 * t304) * t178;
t62 = (-t105 + t240) * t178 + t230;
t188 = (t62 * (-t153 - t168 - t172) - t256 + (t62 * (-rSges(5,3) - qJ(4)) + t63 * t243) * t170) * t178;
t149 = rSges(4,2) * t273;
t142 = rSges(5,2) * t248;
t125 = t178 * t134;
t114 = t144 * t178 + t251;
t104 = -t127 * t178 - t254;
t103 = -t178 * t274 - t255;
t95 = -t178 * (rSges(4,1) * t271 - t149) + t207;
t94 = -t134 * t177 + t208;
t93 = rSges(6,1) * t117 - rSges(6,2) * t118;
t92 = -rSges(6,1) * t115 - rSges(6,2) * t116;
t53 = (t178 * t298 - t100 + t142) * t178 + t203;
t52 = t178 * (-t178 * t253 + t262) + t198;
t48 = -rSges(6,3) * t249 + t284;
t36 = t221 * t257 + qJD(3);
t21 = t160 * t48 - t171 * t312 + t195 * t273 + t198;
t16 = ((-t178 * t79 + t49) * t171 + (-t178 * t77 - t48) * t170) * t257;
t1 = [m(3) * (t104 * (-t143 - t290) + t103 * (t144 + t176) + (-t127 - t251 + t114) * t113) + t193 + (t22 * (t199 - t290) + t34 * (t216 - t251) + t21 * (t176 + t200) + (t34 + t293 + t296) * t35 + t302) * m(6) + (t53 * (t194 - t290) + t62 * (t142 - t229) + t52 * (t176 + t306) + t188 + (t62 + t296 + t300) * t63) * m(5) + (t95 * (t204 - t290) + t96 * (t149 - t251) + t94 * (t176 + t299) + t189 + (t125 - t205 + t96 - t252) * t97) * m(4); t193 + (pkin(2) * t279 * t178 + t22 * t199 + t21 * t200 + (t293 + t295) * t35 + (t216 + t294) * t34 + t302) * m(6) + (-(-t306 * t62 - t256) * t178 + t53 * t194 + t62 * t142 + t52 * t306 + t188 + (t295 + t300) * t63) * m(5) + (t97 * t125 - (-t97 * t289 - t299 * t96) * t178 + t96 * t149 + t95 * t204 + t94 * t299 + t189) * m(4) + (-(-t113 * t144 - t114 * t143) * t178 + t103 * t144 - t104 * t143 - t113 * t127 - t114 * t274) * m(3); m(6) * t16; 0.2e1 * (t278 / 0.2e1 + t21 * t291) * m(6) + 0.2e1 * (t53 * t170 / 0.2e1 + t52 * t291) * m(5); -t10 * t249 / 0.2e1 + (-t181 * t51 + t210) * t219 + (-t181 * t19 + t191) * t233 + (qJD(5) * t192 + t160 * t20) * t272 / 0.2e1 + (-t181 * t50 + t211) * t218 + (-t181 * t20 + t192) * t235 - t181 * (qJD(5) * t190 + t37) / 0.2e1 + t160 * (-t181 * t40 + t190) / 0.2e1 + ((t263 * t117 + t264 * t118 + t137 * t270) * t160 + (t117 * t292 + t202 * t118 + t171 * t213) * t257) * t234 + ((-t263 * t115 + t264 * t116 + t137 * t272) * t160 + (-t115 * t292 + t202 * t116 + t170 * t213) * t257) * t236 - t160 * (-t181 * t137 * t160 + ((-t263 * t182 + t264 * t184) * t160 + ((-t182 * t292 + t202 * t184) * t180 - t220 * t181) * qJD(5)) * t180) / 0.2e1 + (qJD(5) * t191 + t160 * t19 + t178 * t9) * t270 / 0.2e1 + ((-t21 * t79 + t22 * t77 + t34 * t49 - t35 * t48) * t181 + (t16 * t221 + t36 * (-t170 * t48 + t171 * t49 - t79 * t271 - t77 * t273) + t222 * t132 + ((t34 * t178 - t21) * t171 + (t35 * t178 + t22) * t170) * t128) * t180 - (-t34 * t92 + t35 * t93) * t160 - (t36 * (-t170 * t93 + t171 * t92) + t222 * t212) * t257) * m(6);];
tauc = t1(:);
