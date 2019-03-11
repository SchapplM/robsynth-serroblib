% Calculate joint inertia matrix for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:06
% EndTime: 2019-03-09 11:43:16
% DurationCPUTime: 4.34s
% Computational Cost: add. (9234->408), mult. (8563->590), div. (0->0), fcn. (9160->10), ass. (0->202)
t201 = qJ(2) + pkin(10);
t194 = qJ(4) + t201;
t189 = cos(t194);
t208 = cos(qJ(5));
t210 = cos(qJ(1));
t261 = t208 * t210;
t205 = sin(qJ(5));
t207 = sin(qJ(1));
t263 = t207 * t205;
t156 = t189 * t263 + t261;
t262 = t207 * t208;
t264 = t205 * t210;
t157 = t189 * t262 - t264;
t188 = sin(t194);
t269 = t188 * t207;
t82 = Icges(7,5) * t157 + Icges(7,6) * t269 + Icges(7,3) * t156;
t88 = Icges(6,4) * t157 - Icges(6,2) * t156 + Icges(6,6) * t269;
t323 = t82 - t88;
t158 = t189 * t264 - t262;
t159 = t189 * t261 + t263;
t267 = t188 * t210;
t83 = Icges(7,5) * t159 + Icges(7,6) * t267 + Icges(7,3) * t158;
t89 = Icges(6,4) * t159 - Icges(6,2) * t158 + Icges(6,6) * t267;
t322 = t83 - t89;
t84 = Icges(6,5) * t157 - Icges(6,6) * t156 + Icges(6,3) * t269;
t86 = Icges(7,4) * t157 + Icges(7,2) * t269 + Icges(7,6) * t156;
t321 = t84 + t86;
t85 = Icges(6,5) * t159 - Icges(6,6) * t158 + Icges(6,3) * t267;
t87 = Icges(7,4) * t159 + Icges(7,2) * t267 + Icges(7,6) * t158;
t320 = t85 + t87;
t90 = Icges(7,1) * t157 + Icges(7,4) * t269 + Icges(7,5) * t156;
t92 = Icges(6,1) * t157 - Icges(6,4) * t156 + Icges(6,5) * t269;
t319 = t90 + t92;
t91 = Icges(7,1) * t159 + Icges(7,4) * t267 + Icges(7,5) * t158;
t93 = Icges(6,1) * t159 - Icges(6,4) * t158 + Icges(6,5) * t267;
t318 = t91 + t93;
t301 = rSges(7,3) + qJ(6);
t305 = rSges(7,1) + pkin(5);
t317 = -t301 * t156 - t305 * t157;
t316 = t156 * t323 + t319 * t157 + t321 * t269;
t315 = t156 * t322 + t157 * t318 + t269 * t320;
t314 = t158 * t323 + t319 * t159 + t321 * t267;
t313 = t158 * t322 + t159 * t318 + t267 * t320;
t114 = -Icges(7,6) * t189 + (Icges(7,5) * t208 + Icges(7,3) * t205) * t188;
t116 = -Icges(7,2) * t189 + (Icges(7,4) * t208 + Icges(7,6) * t205) * t188;
t118 = -Icges(7,4) * t189 + (Icges(7,1) * t208 + Icges(7,5) * t205) * t188;
t51 = t114 * t156 + t116 * t269 + t118 * t157;
t115 = -Icges(6,3) * t189 + (Icges(6,5) * t208 - Icges(6,6) * t205) * t188;
t117 = -Icges(6,6) * t189 + (Icges(6,4) * t208 - Icges(6,2) * t205) * t188;
t119 = -Icges(6,5) * t189 + (Icges(6,1) * t208 - Icges(6,4) * t205) * t188;
t52 = t115 * t269 - t117 * t156 + t119 * t157;
t312 = -t51 - t52;
t53 = t158 * t114 + t116 * t267 + t159 * t118;
t54 = t115 * t267 - t158 * t117 + t159 * t119;
t311 = -t53 - t54;
t278 = rSges(7,2) * t269 - t317;
t310 = rSges(7,2) * t267 + t301 * t158 + t159 * t305;
t309 = Icges(3,3) + Icges(4,3);
t192 = sin(t201);
t193 = cos(t201);
t206 = sin(qJ(2));
t209 = cos(qJ(2));
t308 = Icges(3,5) * t209 + Icges(4,5) * t193 - Icges(3,6) * t206 - Icges(4,6) * t192;
t202 = t207 ^ 2;
t307 = t312 * t189 + (t316 * t207 + t315 * t210) * t188;
t306 = t311 * t189 + (t314 * t207 + t313 * t210) * t188;
t304 = t207 * pkin(7);
t303 = t315 * t207 - t316 * t210;
t302 = t313 * t207 - t314 * t210;
t300 = -t115 - t116;
t240 = rSges(4,1) * t193 - rSges(4,2) * t192;
t238 = -t157 * rSges(6,1) + t156 * rSges(6,2);
t95 = rSges(6,3) * t269 - t238;
t97 = t159 * rSges(6,1) - t158 * rSges(6,2) + rSges(6,3) * t267;
t299 = t207 * t95 + t210 * t97;
t270 = t188 * t205;
t298 = t114 * t270 + (t118 + t119) * t188 * t208;
t297 = -t308 * t207 + t309 * t210;
t296 = t309 * t207 + t308 * t210;
t229 = Icges(5,5) * t189 - Icges(5,6) * t188;
t128 = -Icges(5,3) * t210 + t207 * t229;
t129 = Icges(5,3) * t207 + t210 * t229;
t203 = t210 ^ 2;
t271 = Icges(5,4) * t189;
t232 = -Icges(5,2) * t188 + t271;
t131 = Icges(5,6) * t207 + t210 * t232;
t272 = Icges(5,4) * t188;
t235 = Icges(5,1) * t189 - t272;
t133 = Icges(5,5) * t207 + t210 * t235;
t227 = -t131 * t188 + t133 * t189;
t130 = -Icges(5,6) * t210 + t207 * t232;
t132 = -Icges(5,5) * t210 + t207 * t235;
t228 = t130 * t188 - t132 * t189;
t295 = -t203 * t128 - (t227 * t207 + (-t129 + t228) * t210) * t207 - t303;
t294 = t278 * t207 + t310 * t210;
t292 = t207 / 0.2e1;
t291 = -t210 / 0.2e1;
t290 = pkin(2) * t206;
t289 = pkin(4) * t189;
t204 = -qJ(3) - pkin(7);
t191 = t209 * pkin(2) + pkin(1);
t288 = -t117 * t270 + t189 * t300 + t298;
t287 = rSges(3,1) * t209;
t285 = rSges(3,2) * t206;
t283 = t210 * rSges(3,3);
t40 = -t189 * t86 + (t205 * t82 + t208 * t90) * t188;
t282 = t40 * t210;
t41 = -t189 * t87 + (t205 * t83 + t208 * t91) * t188;
t281 = t41 * t207;
t42 = -t189 * t84 + (-t205 * t88 + t208 * t92) * t188;
t280 = t42 * t210;
t43 = -t189 * t85 + (-t205 * t89 + t208 * t93) * t188;
t279 = t43 * t207;
t276 = Icges(3,4) * t206;
t275 = Icges(3,4) * t209;
t274 = Icges(4,4) * t192;
t273 = Icges(4,4) * t193;
t266 = t189 * t210;
t200 = -pkin(8) + t204;
t265 = t200 * t210;
t259 = -t189 * rSges(7,2) + (t205 * t301 + t208 * t305) * t188;
t121 = -t189 * rSges(6,3) + (rSges(6,1) * t208 - rSges(6,2) * t205) * t188;
t165 = pkin(4) * t188 - pkin(9) * t189;
t258 = -t121 - t165;
t220 = rSges(5,1) * t266 - rSges(5,2) * t267 + t207 * rSges(5,3);
t239 = rSges(5,1) * t189 - rSges(5,2) * t188;
t75 = t207 * (-t210 * rSges(5,3) + t207 * t239) + t210 * t220;
t185 = t210 * t191;
t199 = t210 * pkin(7);
t257 = t207 * (t199 + (-pkin(1) + t191) * t207) + t210 * (-pkin(1) * t210 + t185 - t304);
t255 = pkin(4) * t266 + pkin(9) * t267;
t256 = t202 * (pkin(9) * t188 + t289) + t210 * t255;
t254 = t207 * rSges(3,3) + t210 * t287;
t252 = t202 + t203;
t251 = (t202 * t129 + (t228 * t210 + (-t128 + t227) * t207) * t210 + t302) * t207;
t250 = -t165 - t259;
t247 = -rSges(4,1) * t192 - rSges(4,2) * t193 - t290;
t171 = pkin(3) * t193 + t191;
t246 = -t171 - t289;
t166 = t210 * t171;
t245 = -t207 * t200 + t166;
t244 = t210 * (t166 - t185) + t257 + (t171 - t191) * t202;
t242 = -pkin(3) * t192 - t290;
t241 = -t285 + t287;
t237 = Icges(3,1) * t209 - t276;
t236 = Icges(4,1) * t193 - t274;
t234 = -Icges(3,2) * t206 + t275;
t233 = -Icges(4,2) * t192 + t273;
t162 = Icges(5,2) * t189 + t272;
t163 = Icges(5,1) * t188 + t271;
t222 = -t162 * t188 + t163 * t189;
t221 = t207 * rSges(4,3) + t210 * t240;
t164 = rSges(5,1) * t188 + rSges(5,2) * t189;
t219 = -t164 + t242;
t218 = -t165 + t242;
t217 = t245 + t255;
t216 = t244 + t256;
t215 = -t121 + t218;
t214 = t295 * t210 + t251;
t213 = t218 - t259;
t212 = -(t281 - t282 + t279 - t280) * t189 / 0.2e1 + t306 * t292 + t307 * t291 + t303 * t269 / 0.2e1 + t302 * t267 / 0.2e1;
t161 = Icges(5,5) * t188 + Icges(5,6) * t189;
t211 = -t282 / 0.2e1 + t281 / 0.2e1 - t280 / 0.2e1 + t279 / 0.2e1 + (t131 * t189 + t133 * t188 + t207 * t161 + t210 * t222 - t311) * t292 + (t130 * t189 + t132 * t188 - t161 * t210 + t207 * t222 - t312) * t291;
t183 = rSges(2,1) * t210 - t207 * rSges(2,2);
t182 = -t207 * rSges(2,1) - rSges(2,2) * t210;
t181 = rSges(3,1) * t206 + rSges(3,2) * t209;
t135 = t247 * t210;
t134 = t247 * t207;
t127 = t304 + (pkin(1) - t285) * t210 + t254;
t126 = t283 + t199 + (-pkin(1) - t241) * t207;
t113 = t219 * t210;
t112 = t219 * t207;
t111 = -t207 * t204 + t185 + t221;
t110 = (rSges(4,3) - t204) * t210 + (-t191 - t240) * t207;
t100 = t210 * (-t210 * t285 + t254) + (t207 * t241 - t283) * t207;
t99 = t220 + t245;
t98 = (rSges(5,3) - t200) * t210 + (-t171 - t239) * t207;
t81 = t258 * t210;
t80 = t258 * t207;
t74 = t215 * t210;
t73 = t215 * t207;
t68 = t250 * t210;
t67 = t250 * t207;
t66 = t213 * t210;
t65 = t213 * t207;
t64 = t217 + t97;
t63 = -t265 + ((-rSges(6,3) - pkin(9)) * t188 + t246) * t207 + t238;
t62 = -t121 * t267 - t189 * t97;
t61 = t121 * t269 + t189 * t95;
t60 = t210 * t221 + (-t210 * rSges(4,3) + t207 * t240) * t207 + t257;
t57 = (-t207 * t97 + t210 * t95) * t188;
t56 = t217 + t310;
t55 = -t265 + ((-rSges(7,2) - pkin(9)) * t188 + t246) * t207 + t317;
t46 = t256 + t299;
t45 = -t189 * t310 - t259 * t267;
t44 = t189 * t278 + t259 * t269;
t35 = t244 + t75;
t24 = (-t207 * t310 + t210 * t278) * t188;
t23 = t256 + t294;
t22 = t216 + t299;
t19 = t216 + t294;
t1 = [t193 * (Icges(4,2) * t193 + t274) + t192 * (Icges(4,1) * t192 + t273) + t209 * (Icges(3,2) * t209 + t276) + t206 * (Icges(3,1) * t206 + t275) + Icges(2,3) + (-t117 * t205 + t163) * t188 + (t162 + t300) * t189 + m(7) * (t55 ^ 2 + t56 ^ 2) + m(6) * (t63 ^ 2 + t64 ^ 2) + m(5) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + m(4) * (t110 ^ 2 + t111 ^ 2) + m(2) * (t182 ^ 2 + t183 ^ 2) + t298; m(3) * (-t126 * t210 - t127 * t207) * t181 + m(7) * (t55 * t66 + t56 * t65) + m(6) * (t63 * t74 + t64 * t73) + m(5) * (t112 * t99 + t113 * t98) + m(4) * (t110 * t135 + t111 * t134) + t211 + (t193 * (Icges(4,6) * t207 + t210 * t233) + t192 * (Icges(4,5) * t207 + t210 * t236) + t209 * (Icges(3,6) * t207 + t210 * t234) + t206 * (Icges(3,5) * t207 + t210 * t237)) * t292 + (t193 * (-Icges(4,6) * t210 + t207 * t233) + t192 * (-Icges(4,5) * t210 + t207 * t236) + t209 * (-Icges(3,6) * t210 + t207 * t234) + t206 * (-Icges(3,5) * t210 + t207 * t237)) * t291 + (Icges(3,5) * t206 + Icges(4,5) * t192 + Icges(3,6) * t209 + Icges(4,6) * t193) * (t202 / 0.2e1 + t203 / 0.2e1); m(7) * (t19 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t22 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t35 ^ 2) + m(4) * (t134 ^ 2 + t135 ^ 2 + t60 ^ 2) + m(3) * (t181 ^ 2 * t252 + t100 ^ 2) + t251 + t296 * t207 * t202 + (t297 * t203 + (t297 * t207 + t296 * t210) * t207 + t295) * t210; m(7) * (t207 * t55 - t210 * t56) + m(6) * (t207 * t63 - t210 * t64) + m(5) * (t207 * t98 - t210 * t99) + m(4) * (t207 * t110 - t111 * t210); m(7) * (t207 * t66 - t210 * t65) + m(6) * (t207 * t74 - t210 * t73) + m(5) * (-t112 * t210 + t207 * t113) + m(4) * (-t134 * t210 + t207 * t135); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t252; m(7) * (t55 * t68 + t56 * t67) + m(6) * (t63 * t81 + t64 * t80) + m(5) * (-t207 * t99 - t210 * t98) * t164 + t211; m(7) * (t19 * t23 + t65 * t67 + t66 * t68) + m(6) * (t22 * t46 + t73 * t80 + t74 * t81) + m(5) * (t75 * t35 + (-t112 * t207 - t113 * t210) * t164) + t214; m(6) * (t81 * t207 - t210 * t80) + m(7) * (t68 * t207 - t210 * t67); m(7) * (t23 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t46 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t164 ^ 2 * t252 + t75 ^ 2) + t214; -t288 * t189 + m(7) * (t44 * t55 + t45 * t56) + m(6) * (t61 * t63 + t62 * t64) + ((t43 / 0.2e1 + t54 / 0.2e1 + t53 / 0.2e1 + t41 / 0.2e1) * t210 + (t42 / 0.2e1 + t52 / 0.2e1 + t51 / 0.2e1 + t40 / 0.2e1) * t207) * t188; m(7) * (t19 * t24 + t44 * t66 + t45 * t65) + m(6) * (t22 * t57 + t61 * t74 + t62 * t73) + t212; m(6) * (t61 * t207 - t210 * t62) + m(7) * (t44 * t207 - t210 * t45); m(7) * (t23 * t24 + t44 * t68 + t45 * t67) + m(6) * (t46 * t57 + t61 * t81 + t62 * t80) + t212; m(7) * (t24 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t57 ^ 2 + t61 ^ 2 + t62 ^ 2) + t288 * t189 ^ 2 + (t306 * t210 + t307 * t207 + ((-t41 - t43) * t210 + (-t40 - t42) * t207) * t189) * t188; m(7) * (t156 * t56 + t158 * t55); m(7) * (t156 * t65 + t158 * t66 + t19 * t270); m(7) * (-t156 * t210 + t158 * t207); m(7) * (t156 * t67 + t158 * t68 + t23 * t270); m(7) * (t156 * t45 + t158 * t44 + t24 * t270); m(7) * (t188 ^ 2 * t205 ^ 2 + t156 ^ 2 + t158 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
