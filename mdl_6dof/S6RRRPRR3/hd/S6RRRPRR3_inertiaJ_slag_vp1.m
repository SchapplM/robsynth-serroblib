% Calculate joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:11:06
% EndTime: 2019-03-09 18:11:16
% DurationCPUTime: 4.92s
% Computational Cost: add. (11317->441), mult. (14895->636), div. (0->0), fcn. (17358->10), ass. (0->221)
t217 = cos(qJ(1));
t212 = sin(qJ(5));
t210 = qJ(2) + qJ(3);
t202 = sin(t210);
t322 = cos(qJ(5));
t269 = t202 * t322;
t203 = cos(t210);
t283 = t203 * t217;
t161 = t212 * t283 - t217 * t269;
t171 = t202 * t212 + t203 * t322;
t162 = t171 * t217;
t214 = sin(qJ(1));
t102 = Icges(6,5) * t162 - Icges(6,6) * t161 - Icges(6,3) * t214;
t103 = Icges(6,4) * t162 - Icges(6,2) * t161 - Icges(6,6) * t214;
t104 = Icges(6,1) * t162 - Icges(6,4) * t161 - Icges(6,5) * t214;
t160 = t171 * t214;
t211 = sin(qJ(6));
t215 = cos(qJ(6));
t126 = -t160 * t211 + t215 * t217;
t127 = t160 * t215 + t211 * t217;
t284 = t203 * t214;
t159 = t212 * t284 - t214 * t269;
t287 = Icges(7,6) * t159;
t289 = Icges(7,2) * t126;
t291 = Icges(7,5) * t159;
t295 = Icges(7,4) * t127;
t301 = Icges(7,1) * t127;
t220 = Icges(7,3) * t159 ^ 2 + (0.2e1 * t291 + t301) * t127 + (0.2e1 * t287 + t289 + 0.2e1 * t295) * t126;
t128 = -t162 * t211 - t214 * t215;
t129 = t162 * t215 - t211 * t214;
t74 = Icges(7,5) * t129 + Icges(7,6) * t128 + Icges(7,3) * t161;
t75 = Icges(7,4) * t129 + Icges(7,2) * t128 + Icges(7,6) * t161;
t76 = Icges(7,1) * t129 + Icges(7,4) * t128 + Icges(7,5) * t161;
t27 = t126 * t75 + t127 * t76 + t159 * t74;
t12 = t27 * t214 - t217 * t220;
t209 = t217 ^ 2;
t288 = Icges(6,6) * t217;
t290 = Icges(6,2) * t159;
t292 = Icges(6,5) * t217;
t296 = Icges(6,4) * t160;
t302 = Icges(6,1) * t160;
t336 = (t102 * t217 - t103 * t159 + t104 * t160) * t214 - (Icges(6,3) * t209 + (0.2e1 * t292 + t302) * t160 + (-0.2e1 * t288 + t290 - 0.2e1 * t296) * t159) * t217 + t12;
t329 = t336 * t217;
t224 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t159;
t225 = t287 + t289 + t295;
t227 = Icges(7,4) * t126 + t291 + t301;
t219 = t128 * t225 + t129 * t227 + t161 * t224;
t28 = t128 * t75 + t129 * t76 + t161 * t74;
t13 = t28 * t214 - t217 * t219;
t226 = t288 - t290 + t296;
t228 = -Icges(6,4) * t159 + t292 + t302;
t331 = ((-t102 * t214 - t103 * t161 + t104 * t162) * t214 - (t162 * t228 - t161 * t226 - t214 * (Icges(6,5) * t160 - Icges(6,6) * t159 + Icges(6,3) * t217)) * t217 + t13) * t214;
t231 = t329 - t331;
t294 = Icges(5,5) * t202;
t298 = Icges(4,4) * t202;
t338 = t294 - t298 + (-Icges(4,2) - Icges(5,3)) * t203;
t293 = Icges(5,5) * t203;
t297 = Icges(4,4) * t203;
t337 = -t293 + t297 + (Icges(4,1) + Icges(5,1)) * t202;
t324 = t214 / 0.2e1;
t323 = -t217 / 0.2e1;
t335 = t214 * pkin(7);
t78 = t129 * rSges(7,1) + t128 * rSges(7,2) + t161 * rSges(7,3);
t334 = -t162 * pkin(5) - t161 * pkin(10) - t78;
t333 = (Icges(4,6) - Icges(5,6)) * t203 + (Icges(5,4) + Icges(4,5)) * t202;
t332 = t202 * t338 + t337 * t203;
t245 = Icges(4,5) * t203 - Icges(4,6) * t202;
t142 = -Icges(4,3) * t217 + t245 * t214;
t143 = Icges(4,3) * t214 + t245 * t217;
t247 = Icges(5,4) * t203 + Icges(5,6) * t202;
t144 = -Icges(5,2) * t217 + t247 * t214;
t145 = Icges(5,2) * t214 + t247 * t217;
t248 = -Icges(4,2) * t202 + t297;
t147 = Icges(4,6) * t214 + t248 * t217;
t251 = Icges(4,1) * t203 - t298;
t151 = Icges(4,5) * t214 + t251 * t217;
t240 = -t147 * t202 + t151 * t203;
t146 = -Icges(4,6) * t217 + t248 * t214;
t150 = -Icges(4,5) * t217 + t251 * t214;
t241 = t146 * t202 - t150 * t203;
t244 = Icges(5,3) * t202 + t293;
t141 = Icges(5,6) * t214 + t244 * t217;
t250 = Icges(5,1) * t203 + t294;
t149 = Icges(5,4) * t214 + t250 * t217;
t242 = t141 * t202 + t149 * t203;
t140 = -Icges(5,6) * t217 + t244 * t214;
t148 = -Icges(5,4) * t217 + t250 * t214;
t243 = -t140 * t202 - t148 * t203;
t330 = (-t144 - t142) * t209 + ((-t242 - t240) * t214 + (t145 - t243 + t143 - t241) * t217) * t214;
t172 = -t203 * t212 + t269;
t93 = Icges(7,3) * t171 + (Icges(7,5) * t215 - Icges(7,6) * t211) * t172;
t94 = Icges(7,6) * t171 + (Icges(7,4) * t215 - Icges(7,2) * t211) * t172;
t95 = Icges(7,5) * t171 + (Icges(7,1) * t215 - Icges(7,4) * t211) * t172;
t36 = t126 * t94 + t127 * t95 + t159 * t93;
t3 = t159 * t220 + t27 * t161 + t36 * t171;
t32 = t171 * t74 + (-t211 * t75 + t215 * t76) * t172;
t304 = t32 * t214;
t31 = t171 * t224 + (-t211 * t225 + t215 * t227) * t172;
t305 = t31 * t217;
t37 = t128 * t94 + t129 * t95 + t161 * t93;
t4 = t159 * t219 + t28 * t161 + t37 * t171;
t263 = t4 * t324 + t171 * (t304 - t305) / 0.2e1 + t3 * t323 + t159 * t12 / 0.2e1 + t161 * t13 / 0.2e1;
t208 = t214 ^ 2;
t328 = m(5) / 0.2e1;
t327 = m(7) / 0.2e1;
t213 = sin(qJ(2));
t321 = pkin(2) * t213;
t320 = t160 * pkin(5);
t218 = -pkin(8) - pkin(7);
t317 = -pkin(9) - t218;
t315 = t172 * t215 * t95 + t171 * t93;
t216 = cos(qJ(2));
t314 = rSges(3,1) * t216;
t313 = rSges(3,2) * t213;
t309 = t211 * t94;
t308 = t217 * rSges(3,3);
t96 = t171 * rSges(7,3) + (rSges(7,1) * t215 - rSges(7,2) * t211) * t172;
t303 = pkin(5) * t172 + pkin(10) * t171 + t96;
t256 = -t160 * rSges(6,1) + t159 * rSges(6,2);
t281 = t162 * rSges(6,1) - t161 * rSges(6,2);
t65 = t214 * t256 - t217 * t281;
t300 = Icges(3,4) * t213;
t299 = Icges(3,4) * t216;
t286 = qJ(4) * t202;
t285 = t202 * t217;
t200 = pkin(2) * t216 + pkin(1);
t193 = t217 * t200;
t207 = t217 * pkin(7);
t282 = t214 * (t207 + (-pkin(1) + t200) * t214) + t217 * (-t217 * pkin(1) + t193 - t335);
t234 = rSges(4,1) * t283 - rSges(4,2) * t285 + t214 * rSges(4,3);
t257 = rSges(4,1) * t203 - rSges(4,2) * t202;
t101 = t214 * (-t217 * rSges(4,3) + t257 * t214) + t217 * t234;
t278 = pkin(3) * t283 + qJ(4) * t285;
t280 = t208 * (pkin(3) * t203 + t286) + t217 * t278;
t180 = pkin(3) * t202 - qJ(4) * t203;
t279 = -rSges(5,1) * t202 + rSges(5,3) * t203 - t180;
t277 = t214 * rSges(3,3) + t217 * t314;
t276 = t208 + t209;
t275 = -rSges(6,3) + t317;
t274 = t32 / 0.2e1 + t37 / 0.2e1;
t273 = t36 / 0.2e1 + t31 / 0.2e1;
t270 = rSges(5,1) * t283 + t214 * rSges(5,2) + rSges(5,3) * t285;
t182 = rSges(4,1) * t202 + rSges(4,2) * t203;
t268 = -t182 - t321;
t267 = -pkin(4) * t202 - t180;
t265 = t214 * (t208 * t145 + (t243 * t217 + (-t144 + t242) * t214) * t217) + t214 * (t208 * t143 + (t241 * t217 + (-t142 + t240) * t214) * t217) + t331;
t264 = -t214 * t218 + t193;
t255 = -t127 * rSges(7,1) - t126 * rSges(7,2);
t77 = t159 * rSges(7,3) - t255;
t33 = t334 * t217 + (-t159 * pkin(10) - t320 - t77) * t214;
t71 = t214 * (-t217 * rSges(5,2) + (rSges(5,1) * t203 + rSges(5,3) * t202) * t214) + t217 * t270 + t280;
t197 = pkin(4) * t283;
t262 = t214 * pkin(4) * t284 + t217 * t197 + t280;
t261 = t193 + t197 + t278;
t260 = t279 - t321;
t120 = rSges(6,1) * t172 - rSges(6,2) * t171;
t259 = -t120 + t267;
t258 = -t313 + t314;
t233 = t267 - t321;
t230 = -t120 + t233;
t81 = t230 * t214;
t82 = t230 * t217;
t254 = t214 * t81 + t217 * t82;
t83 = t259 * t214;
t84 = t259 * t217;
t253 = t214 * t83 + t217 * t84;
t252 = Icges(3,1) * t216 - t300;
t249 = -Icges(3,2) * t213 + t299;
t246 = Icges(3,5) * t216 - Icges(3,6) * t213;
t235 = t267 - t303;
t47 = t262 - t65;
t223 = (-t286 - t200 + (-pkin(3) - pkin(4)) * t203) * t214;
t69 = t275 * t217 + t223 + t256;
t70 = t275 * t214 + t261 + t281;
t232 = m(6) * (t214 * t70 + t217 * t69);
t229 = t233 - t303;
t24 = t262 - t33;
t222 = (-t336 + t330) * t217 + t265;
t117 = Icges(6,5) * t172 - Icges(6,6) * t171;
t118 = Icges(6,4) * t172 - Icges(6,2) * t171;
t119 = Icges(6,1) * t172 - Icges(6,4) * t171;
t54 = t117 * t217 - t118 * t159 + t119 * t160;
t55 = -t117 * t214 - t118 * t161 + t119 * t162;
t59 = -t171 * t226 + t172 * t228;
t60 = -t103 * t171 + t104 * t172;
t221 = -t305 / 0.2e1 + t304 / 0.2e1 + (t37 + t55 + t60 + t332 * t217 + t333 * t214 + (-t141 + t147) * t203 + (t149 + t151) * t202) * t324 + (t36 + t54 + t59 - t333 * t217 + t332 * t214 + (-t140 + t146) * t203 + (t148 + t150) * t202) * t323;
t191 = rSges(2,1) * t217 - rSges(2,2) * t214;
t190 = -rSges(2,1) * t214 - rSges(2,2) * t217;
t189 = rSges(3,1) * t213 + rSges(3,2) * t216;
t164 = Icges(3,3) * t214 + t246 * t217;
t163 = -Icges(3,3) * t217 + t246 * t214;
t139 = t268 * t217;
t138 = t268 * t214;
t131 = t335 + (pkin(1) - t313) * t217 + t277;
t130 = t308 + t207 + (-pkin(1) - t258) * t214;
t123 = t279 * t217;
t122 = t279 * t214;
t115 = t234 + t264;
t114 = (rSges(4,3) - t218) * t217 + (-t200 - t257) * t214;
t113 = t260 * t217;
t112 = t260 * t214;
t111 = t217 * (-t217 * t313 + t277) + (t258 * t214 - t308) * t214;
t98 = t264 + t270 + t278;
t97 = (rSges(5,2) - t218) * t217 + (-t200 + (-rSges(5,1) - pkin(3)) * t203 + (-rSges(5,3) - qJ(4)) * t202) * t214;
t68 = t303 * t217;
t67 = t303 * t214;
t66 = t101 + t282;
t64 = t235 * t217;
t63 = t235 * t214;
t62 = t229 * t217;
t61 = t229 * t214;
t56 = t71 + t282;
t51 = t317 * t214 + t261 - t334;
t50 = -t320 + t317 * t217 + (-rSges(7,3) - pkin(10)) * t159 + t223 + t255;
t44 = -t161 * t96 + t171 * t78;
t43 = t159 * t96 - t171 * t77;
t40 = -t159 * t78 + t161 * t77;
t39 = t47 + t282;
t38 = (-t172 * t309 + t315) * t171;
t21 = t24 + t282;
t1 = [-t171 * t118 + t216 * (Icges(3,2) * t216 + t300) + t213 * (Icges(3,1) * t213 + t299) + Icges(2,3) - t338 * t203 + t337 * t202 + (t119 - t309) * t172 + m(7) * (t50 ^ 2 + t51 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t97 ^ 2 + t98 ^ 2) + m(4) * (t114 ^ 2 + t115 ^ 2) + m(3) * (t130 ^ 2 + t131 ^ 2) + m(2) * (t190 ^ 2 + t191 ^ 2) + t315; (t208 / 0.2e1 + t209 / 0.2e1) * (Icges(3,5) * t213 + Icges(3,6) * t216) + m(7) * (t50 * t62 + t51 * t61) + m(6) * (t69 * t82 + t70 * t81) + m(5) * (t112 * t98 + t113 * t97) + m(4) * (t114 * t139 + t115 * t138) + t221 + m(3) * (-t130 * t217 - t131 * t214) * t189 + (t216 * (Icges(3,6) * t214 + t249 * t217) + t213 * (Icges(3,5) * t214 + t252 * t217)) * t324 + (t216 * (-Icges(3,6) * t217 + t249 * t214) + t213 * (-Icges(3,5) * t217 + t252 * t214)) * t323; m(7) * (t21 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t39 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t56 ^ 2) + m(4) * (t138 ^ 2 + t139 ^ 2 + t66 ^ 2) + m(3) * (t189 ^ 2 * t276 + t111 ^ 2) + t214 * t208 * t164 + t265 + (-t209 * t163 + (-t214 * t163 + t217 * t164) * t214 + t330) * t217 - t329; m(7) * (t50 * t64 + t51 * t63) + m(6) * (t69 * t84 + t70 * t83) + m(5) * (t122 * t98 + t123 * t97) + m(4) * (-t114 * t217 - t115 * t214) * t182 + t221; m(7) * (t21 * t24 + t61 * t63 + t62 * t64) + m(6) * (t39 * t47 + t81 * t83 + t82 * t84) + m(5) * (t112 * t122 + t113 * t123 + t56 * t71) + m(4) * (t101 * t66 + (-t138 * t214 - t139 * t217) * t182) + t222; m(7) * (t24 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t47 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(5) * (t122 ^ 2 + t123 ^ 2 + t71 ^ 2) + m(4) * (t182 ^ 2 * t276 + t101 ^ 2) + t222; 0.2e1 * ((t214 * t51 + t217 * t50) * t327 + t232 / 0.2e1 + (t214 * t98 + t217 * t97) * t328) * t202; m(7) * (-t203 * t21 + (t214 * t61 + t217 * t62) * t202) + m(6) * (t202 * t254 - t203 * t39) + m(5) * (-t203 * t56 + (t112 * t214 + t113 * t217) * t202); m(7) * (-t203 * t24 + (t214 * t63 + t217 * t64) * t202) + m(6) * (t202 * t253 - t203 * t47) + m(5) * (-t203 * t71 + (t122 * t214 + t123 * t217) * t202); 0.2e1 * (t328 + m(6) / 0.2e1 + t327) * (t202 ^ 2 * t276 + t203 ^ 2); (t59 / 0.2e1 + t54 / 0.2e1 + t273) * t217 + (-t55 / 0.2e1 - t60 / 0.2e1 - t274) * t214 + m(7) * (t50 * t68 + t51 * t67) + t120 * t232; m(7) * (t21 * t33 + t61 * t67 + t62 * t68) + m(6) * (t120 * t254 + t65 * t39) + t231; m(7) * (t24 * t33 + t63 * t67 + t64 * t68) + m(6) * (t120 * t253 + t65 * t47) + t231; m(6) * (t120 * t202 * t276 - t65 * t203) + m(7) * (-t33 * t203 + (t214 * t67 + t217 * t68) * t202); m(6) * (t120 ^ 2 * t276 + t65 ^ 2) + m(7) * (t33 ^ 2 + t67 ^ 2 + t68 ^ 2) - t231; m(7) * (t43 * t50 + t44 * t51) + t38 + t274 * t161 + t273 * t159; m(7) * (t21 * t40 + t43 * t62 + t44 * t61) + t263; m(7) * (t24 * t40 + t43 * t64 + t44 * t63) + t263; m(7) * (-t40 * t203 + (t214 * t44 + t217 * t43) * t202); m(7) * (t33 * t40 + t43 * t68 + t44 * t67) - t263; t161 * t4 + t159 * t3 + t171 * (t31 * t159 + t32 * t161 + t38) + m(7) * (t40 ^ 2 + t43 ^ 2 + t44 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
