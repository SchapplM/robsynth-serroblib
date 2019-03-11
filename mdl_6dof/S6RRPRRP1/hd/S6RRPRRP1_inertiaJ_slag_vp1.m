% Calculate joint inertia matrix for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:01
% EndTime: 2019-03-09 11:39:11
% DurationCPUTime: 4.53s
% Computational Cost: add. (9697->417), mult. (8612->596), div. (0->0), fcn. (9102->10), ass. (0->204)
t202 = qJ(2) + pkin(10);
t195 = qJ(4) + t202;
t190 = cos(t195);
t210 = cos(qJ(5));
t212 = cos(qJ(1));
t262 = t210 * t212;
t207 = sin(qJ(5));
t209 = sin(qJ(1));
t264 = t209 * t207;
t154 = -t190 * t264 - t262;
t263 = t209 * t210;
t265 = t207 * t212;
t155 = t190 * t263 - t265;
t189 = sin(t195);
t269 = t189 * t209;
t86 = Icges(7,5) * t155 + Icges(7,6) * t154 + Icges(7,3) * t269;
t88 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t269;
t326 = t86 + t88;
t156 = -t190 * t265 + t263;
t157 = t190 * t262 + t264;
t267 = t189 * t212;
t87 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t267;
t89 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t267;
t325 = t87 + t89;
t90 = Icges(7,4) * t155 + Icges(7,2) * t154 + Icges(7,6) * t269;
t92 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t269;
t324 = t90 + t92;
t91 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t267;
t93 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t267;
t323 = t91 + t93;
t94 = Icges(7,1) * t155 + Icges(7,4) * t154 + Icges(7,5) * t269;
t96 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t269;
t322 = t94 + t96;
t95 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t267;
t97 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t267;
t321 = t95 + t97;
t205 = -qJ(6) - pkin(9);
t320 = rSges(7,3) - t205;
t319 = t324 * t154 + t322 * t155 + t326 * t269;
t318 = t323 * t154 + t321 * t155 + t325 * t269;
t317 = t324 * t156 + t322 * t157 + t326 * t267;
t316 = t323 * t156 + t321 * t157 + t325 * t267;
t114 = -Icges(7,3) * t190 + (Icges(7,5) * t210 - Icges(7,6) * t207) * t189;
t116 = -Icges(7,6) * t190 + (Icges(7,4) * t210 - Icges(7,2) * t207) * t189;
t118 = -Icges(7,5) * t190 + (Icges(7,1) * t210 - Icges(7,4) * t207) * t189;
t51 = t114 * t269 + t116 * t154 + t118 * t155;
t115 = -Icges(6,3) * t190 + (Icges(6,5) * t210 - Icges(6,6) * t207) * t189;
t117 = -Icges(6,6) * t190 + (Icges(6,4) * t210 - Icges(6,2) * t207) * t189;
t119 = -Icges(6,5) * t190 + (Icges(6,1) * t210 - Icges(6,4) * t207) * t189;
t52 = t115 * t269 + t117 * t154 + t119 * t155;
t315 = -t51 - t52;
t53 = t114 * t267 + t156 * t116 + t157 * t118;
t54 = t115 * t267 + t156 * t117 + t157 * t119;
t314 = -t53 - t54;
t191 = pkin(5) * t210 + pkin(4);
t266 = t190 * t212;
t313 = t157 * rSges(7,1) + t156 * rSges(7,2) + pkin(5) * t264 + t191 * t266 + t320 * t267;
t240 = -t155 * rSges(7,1) - t154 * rSges(7,2);
t289 = pkin(9) + t205;
t290 = -pkin(4) + t191;
t287 = -pkin(5) * t265 + (-t189 * t289 + t190 * t290) * t209 + rSges(7,3) * t269 - t240;
t257 = pkin(4) * t266 + pkin(9) * t267;
t312 = -t257 + t313;
t311 = Icges(3,3) + Icges(4,3);
t193 = sin(t202);
t194 = cos(t202);
t208 = sin(qJ(2));
t211 = cos(qJ(2));
t310 = Icges(3,5) * t211 + Icges(4,5) * t194 - Icges(3,6) * t208 - Icges(4,6) * t193;
t203 = t209 ^ 2;
t309 = t315 * t190 + (t209 * t319 + t318 * t212) * t189;
t308 = t314 * t190 + (t209 * t317 + t212 * t316) * t189;
t307 = t209 * pkin(7);
t306 = t318 * t209 - t212 * t319;
t305 = t209 * t316 - t212 * t317;
t304 = -t114 - t115;
t243 = rSges(4,1) * t194 - rSges(4,2) * t193;
t303 = (-t116 - t117) * t207;
t101 = t157 * rSges(6,1) + t156 * rSges(6,2) + rSges(6,3) * t267;
t241 = -t155 * rSges(6,1) - t154 * rSges(6,2);
t99 = rSges(6,3) * t269 - t241;
t302 = t212 * t101 + t209 * t99;
t301 = (t118 + t119) * t189 * t210;
t300 = t311 * t209 + t310 * t212;
t299 = -t310 * t209 + t311 * t212;
t231 = Icges(5,5) * t190 - Icges(5,6) * t189;
t128 = -Icges(5,3) * t212 + t209 * t231;
t129 = Icges(5,3) * t209 + t212 * t231;
t204 = t212 ^ 2;
t271 = Icges(5,4) * t190;
t234 = -Icges(5,2) * t189 + t271;
t131 = Icges(5,6) * t209 + t212 * t234;
t272 = Icges(5,4) * t189;
t237 = Icges(5,1) * t190 - t272;
t133 = Icges(5,5) * t209 + t212 * t237;
t229 = -t131 * t189 + t133 * t190;
t130 = -Icges(5,6) * t212 + t209 * t234;
t132 = -Icges(5,5) * t212 + t209 * t237;
t230 = t130 * t189 - t132 * t190;
t298 = -t204 * t128 - (t229 * t209 + (-t129 + t230) * t212) * t209 - t306;
t297 = t287 * t209 + t312 * t212;
t296 = t190 ^ 2;
t294 = t209 / 0.2e1;
t293 = -t212 / 0.2e1;
t292 = pkin(2) * t208;
t291 = pkin(4) * t190;
t206 = -qJ(3) - pkin(7);
t192 = t211 * pkin(2) + pkin(1);
t288 = t189 * t303 + t190 * t304 + t301;
t286 = rSges(3,1) * t211;
t284 = rSges(3,2) * t208;
t282 = t212 * rSges(3,3);
t42 = -t190 * t86 + (-t207 * t90 + t210 * t94) * t189;
t281 = t42 * t212;
t43 = -t190 * t87 + (-t207 * t91 + t210 * t95) * t189;
t280 = t43 * t209;
t44 = -t190 * t88 + (-t207 * t92 + t210 * t96) * t189;
t279 = t44 * t212;
t45 = -t190 * t89 + (-t207 * t93 + t210 * t97) * t189;
t278 = t45 * t209;
t276 = Icges(3,4) * t208;
t275 = Icges(3,4) * t211;
t274 = Icges(4,4) * t193;
t273 = Icges(4,4) * t194;
t261 = (t289 - rSges(7,3)) * t190 + (rSges(7,1) * t210 - rSges(7,2) * t207 + t290) * t189;
t121 = -rSges(6,3) * t190 + (rSges(6,1) * t210 - rSges(6,2) * t207) * t189;
t163 = pkin(4) * t189 - pkin(9) * t190;
t260 = -t121 - t163;
t222 = rSges(5,1) * t266 - rSges(5,2) * t267 + t209 * rSges(5,3);
t242 = rSges(5,1) * t190 - rSges(5,2) * t189;
t75 = t209 * (-rSges(5,3) * t212 + t209 * t242) + t212 * t222;
t184 = t212 * t192;
t200 = t212 * pkin(7);
t259 = t209 * (t200 + (-pkin(1) + t192) * t209) + t212 * (-pkin(1) * t212 + t184 - t307);
t258 = t203 * (pkin(9) * t189 + t291) + t212 * t257;
t256 = t209 * rSges(3,3) + t212 * t286;
t254 = t203 + t204;
t253 = (t203 * t129 + (t230 * t212 + (-t128 + t229) * t209) * t212 + t305) * t209;
t252 = -t163 - t261;
t249 = -rSges(4,1) * t193 - rSges(4,2) * t194 - t292;
t170 = pkin(3) * t194 + t192;
t164 = t212 * t170;
t201 = -pkin(8) + t206;
t248 = -t209 * t201 + t164;
t247 = t212 * (t164 - t184) + t259 + (t170 - t192) * t203;
t245 = -pkin(3) * t193 - t292;
t244 = -t284 + t286;
t239 = Icges(3,1) * t211 - t276;
t238 = Icges(4,1) * t194 - t274;
t236 = -Icges(3,2) * t208 + t275;
t235 = -Icges(4,2) * t193 + t273;
t160 = Icges(5,2) * t190 + t272;
t161 = Icges(5,1) * t189 + t271;
t224 = -t160 * t189 + t161 * t190;
t223 = t209 * rSges(4,3) + t212 * t243;
t162 = rSges(5,1) * t189 + rSges(5,2) * t190;
t221 = -t162 + t245;
t220 = -t163 + t245;
t218 = t247 + t258;
t217 = -t121 + t220;
t216 = t298 * t212 + t253;
t215 = t220 - t261;
t214 = -(t280 - t281 + t278 - t279) * t190 / 0.2e1 + t308 * t294 + t309 * t293 + t306 * t269 / 0.2e1 + t305 * t267 / 0.2e1;
t159 = Icges(5,5) * t189 + Icges(5,6) * t190;
t213 = -t281 / 0.2e1 + t280 / 0.2e1 - t279 / 0.2e1 + t278 / 0.2e1 + (t131 * t190 + t133 * t189 + t209 * t159 + t212 * t224 - t314) * t294 + (t130 * t190 + t132 * t189 - t159 * t212 + t209 * t224 - t315) * t293;
t182 = rSges(2,1) * t212 - t209 * rSges(2,2);
t181 = -t209 * rSges(2,1) - rSges(2,2) * t212;
t180 = rSges(3,1) * t208 + rSges(3,2) * t211;
t135 = t249 * t212;
t134 = t249 * t209;
t127 = t307 + (pkin(1) - t284) * t212 + t256;
t126 = t282 + t200 + (-pkin(1) - t244) * t209;
t113 = t221 * t212;
t112 = t221 * t209;
t110 = -t209 * t206 + t184 + t223;
t109 = (rSges(4,3) - t206) * t212 + (-t192 - t243) * t209;
t104 = t212 * (-t212 * t284 + t256) + (t209 * t244 - t282) * t209;
t103 = t222 + t248;
t102 = (rSges(5,3) - t201) * t212 + (-t170 - t242) * t209;
t85 = t260 * t212;
t84 = t260 * t209;
t74 = t217 * t212;
t73 = t217 * t209;
t68 = t252 * t212;
t67 = t252 * t209;
t66 = t248 + t101 + t257;
t65 = -t201 * t212 + (-t291 - t170 + (-rSges(6,3) - pkin(9)) * t189) * t209 + t241;
t64 = -t190 * t101 - t121 * t267;
t63 = t121 * t269 + t190 * t99;
t62 = t248 + t313;
t61 = (pkin(5) * t207 - t201) * t212 + (-t320 * t189 - t190 * t191 - t170) * t209 + t240;
t60 = t212 * t223 + (-t212 * rSges(4,3) + t209 * t243) * t209 + t259;
t59 = t215 * t212;
t58 = t215 * t209;
t55 = (-t101 * t209 + t212 * t99) * t189;
t46 = t258 + t302;
t37 = t247 + t75;
t28 = -t190 * t312 - t261 * t267;
t27 = t190 * t287 + t261 * t269;
t24 = (-t209 * t312 + t212 * t287) * t189;
t23 = t258 + t297;
t22 = t218 + t302;
t19 = t218 + t297;
t1 = [t194 * (Icges(4,2) * t194 + t274) + t193 * (Icges(4,1) * t193 + t273) + t211 * (Icges(3,2) * t211 + t276) + t208 * (Icges(3,1) * t208 + t275) + Icges(2,3) + (t160 + t304) * t190 + (t161 + t303) * t189 + m(7) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t65 ^ 2 + t66 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2) + m(4) * (t109 ^ 2 + t110 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + m(2) * (t181 ^ 2 + t182 ^ 2) + t301; m(6) * (t65 * t74 + t66 * t73) + m(7) * (t58 * t62 + t59 * t61) + m(5) * (t102 * t113 + t103 * t112) + m(4) * (t109 * t135 + t110 * t134) + m(3) * (-t126 * t212 - t127 * t209) * t180 + t213 + (t194 * (Icges(4,6) * t209 + t212 * t235) + t193 * (Icges(4,5) * t209 + t212 * t238) + t211 * (Icges(3,6) * t209 + t212 * t236) + t208 * (Icges(3,5) * t209 + t212 * t239)) * t294 + (t194 * (-Icges(4,6) * t212 + t209 * t235) + t193 * (-Icges(4,5) * t212 + t209 * t238) + t211 * (-Icges(3,6) * t212 + t209 * t236) + t208 * (-Icges(3,5) * t212 + t209 * t239)) * t293 + (Icges(3,5) * t208 + Icges(4,5) * t193 + Icges(3,6) * t211 + Icges(4,6) * t194) * (t204 / 0.2e1 + t203 / 0.2e1); m(7) * (t19 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(6) * (t22 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t37 ^ 2) + m(4) * (t134 ^ 2 + t135 ^ 2 + t60 ^ 2) + m(3) * (t180 ^ 2 * t254 + t104 ^ 2) + t253 + t300 * t209 * t203 + (t299 * t204 + (t299 * t209 + t300 * t212) * t209 + t298) * t212; m(7) * (t209 * t61 - t212 * t62) + m(6) * (t209 * t65 - t212 * t66) + m(5) * (t209 * t102 - t103 * t212) + m(4) * (t209 * t109 - t110 * t212); m(7) * (t209 * t59 - t212 * t58) + m(6) * (t209 * t74 - t212 * t73) + m(5) * (-t112 * t212 + t209 * t113) + m(4) * (-t134 * t212 + t209 * t135); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t254; m(7) * (t61 * t68 + t62 * t67) + m(6) * (t65 * t85 + t66 * t84) + m(5) * (-t102 * t212 - t103 * t209) * t162 + t213; m(7) * (t19 * t23 + t58 * t67 + t59 * t68) + m(6) * (t22 * t46 + t73 * t84 + t74 * t85) + m(5) * (t75 * t37 + (-t112 * t209 - t113 * t212) * t162) + t216; m(6) * (t85 * t209 - t212 * t84) + m(7) * (t68 * t209 - t212 * t67); m(7) * (t23 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t46 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t162 ^ 2 * t254 + t75 ^ 2) + t216; -t288 * t190 + m(7) * (t27 * t61 + t28 * t62) + m(6) * (t63 * t65 + t64 * t66) + ((t54 / 0.2e1 + t53 / 0.2e1 + t45 / 0.2e1 + t43 / 0.2e1) * t212 + (t44 / 0.2e1 + t52 / 0.2e1 + t51 / 0.2e1 + t42 / 0.2e1) * t209) * t189; m(7) * (t19 * t24 + t27 * t59 + t28 * t58) + m(6) * (t22 * t55 + t63 * t74 + t64 * t73) + t214; m(6) * (t63 * t209 - t212 * t64) + m(7) * (t27 * t209 - t212 * t28); m(7) * (t23 * t24 + t27 * t68 + t28 * t67) + m(6) * (t46 * t55 + t63 * t85 + t64 * t84) + t214; m(7) * (t24 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t55 ^ 2 + t63 ^ 2 + t64 ^ 2) + t288 * t296 + (t308 * t212 + t309 * t209 + ((-t43 - t45) * t212 + (-t42 - t44) * t209) * t190) * t189; m(7) * (t209 * t62 + t212 * t61) * t189; m(7) * (-t190 * t19 + (t209 * t58 + t212 * t59) * t189); 0; m(7) * (-t190 * t23 + (t209 * t67 + t212 * t68) * t189); m(7) * (-t190 * t24 + (t209 * t28 + t212 * t27) * t189); m(7) * (t189 ^ 2 * t254 + t296);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
