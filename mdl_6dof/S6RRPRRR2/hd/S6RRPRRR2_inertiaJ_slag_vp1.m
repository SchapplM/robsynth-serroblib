% Calculate joint inertia matrix for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:42
% EndTime: 2019-03-09 13:16:52
% DurationCPUTime: 4.23s
% Computational Cost: add. (13803->457), mult. (10647->650), div. (0->0), fcn. (11419->12), ass. (0->232)
t228 = cos(qJ(5));
t207 = pkin(5) * t228 + pkin(4);
t231 = -pkin(10) - pkin(9);
t225 = sin(qJ(5));
t227 = sin(qJ(1));
t287 = t225 * t227;
t220 = qJ(2) + pkin(11);
t211 = qJ(4) + t220;
t206 = cos(t211);
t230 = cos(qJ(1));
t294 = t206 * t230;
t205 = sin(t211);
t295 = t205 * t230;
t223 = qJ(5) + qJ(6);
t213 = cos(t223);
t290 = t213 * t227;
t212 = sin(t223);
t291 = t212 * t230;
t161 = -t206 * t291 + t290;
t289 = t213 * t230;
t292 = t212 * t227;
t162 = t206 * t289 + t292;
t96 = t162 * rSges(7,1) + t161 * rSges(7,2) + rSges(7,3) * t295;
t331 = pkin(5) * t287 + t207 * t294 - t231 * t295 + t96;
t279 = pkin(4) * t294 + pkin(9) * t295;
t330 = -t279 + t331;
t329 = Icges(3,3) + Icges(4,3);
t209 = sin(t220);
t210 = cos(t220);
t226 = sin(qJ(2));
t229 = cos(qJ(2));
t328 = Icges(3,5) * t229 + Icges(4,5) * t210 - Icges(3,6) * t226 - Icges(4,6) * t209;
t221 = t227 ^ 2;
t327 = pkin(7) * t227;
t262 = rSges(4,1) * t210 - rSges(4,2) * t209;
t284 = t228 * t230;
t171 = -t206 * t287 - t284;
t285 = t227 * t228;
t286 = t225 * t230;
t172 = t206 * t285 - t286;
t260 = -rSges(6,1) * t172 - rSges(6,2) * t171;
t296 = t205 * t227;
t112 = rSges(6,3) * t296 - t260;
t173 = -t206 * t286 + t285;
t174 = t206 * t284 + t287;
t113 = t174 * rSges(6,1) + t173 * rSges(6,2) + rSges(6,3) * t295;
t326 = t227 * t112 + t230 * t113;
t325 = t227 * t329 + t230 * t328;
t324 = -t227 * t328 + t230 * t329;
t250 = Icges(5,5) * t206 - Icges(5,6) * t205;
t141 = -Icges(5,3) * t230 + t227 * t250;
t142 = Icges(5,3) * t227 + t230 * t250;
t159 = -t206 * t292 - t289;
t160 = t206 * t290 - t291;
t89 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t296;
t91 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t296;
t93 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t296;
t28 = t159 * t91 + t160 * t93 + t296 * t89;
t90 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t295;
t92 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t295;
t94 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t295;
t29 = t159 * t92 + t160 * t94 + t296 * t90;
t15 = t227 * t29 - t230 * t28;
t106 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t296;
t108 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t296;
t110 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t296;
t42 = t106 * t296 + t108 * t171 + t110 * t172;
t107 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t295;
t109 = Icges(6,4) * t174 + Icges(6,2) * t173 + Icges(6,6) * t295;
t111 = Icges(6,1) * t174 + Icges(6,4) * t173 + Icges(6,5) * t295;
t43 = t107 * t296 + t109 * t171 + t111 * t172;
t22 = t227 * t43 - t230 * t42;
t222 = t230 ^ 2;
t297 = Icges(5,4) * t206;
t253 = -Icges(5,2) * t205 + t297;
t144 = Icges(5,6) * t227 + t230 * t253;
t298 = Icges(5,4) * t205;
t256 = Icges(5,1) * t206 - t298;
t146 = Icges(5,5) * t227 + t230 * t256;
t248 = -t144 * t205 + t146 * t206;
t143 = -Icges(5,6) * t230 + t227 * t253;
t145 = -Icges(5,5) * t230 + t227 * t256;
t249 = t143 * t205 - t145 * t206;
t323 = -t15 - t22 - t222 * t141 - (t248 * t227 + (-t142 + t249) * t230) * t227;
t314 = pkin(9) + t231;
t315 = -pkin(4) + t207;
t102 = -pkin(5) * t286 + (-t205 * t314 + t206 * t315) * t227;
t259 = -t160 * rSges(7,1) - t159 * rSges(7,2);
t95 = rSges(7,3) * t296 - t259;
t322 = t330 * t230 + (t102 + t95) * t227;
t127 = -Icges(7,3) * t206 + (Icges(7,5) * t213 - Icges(7,6) * t212) * t205;
t128 = -Icges(7,6) * t206 + (Icges(7,4) * t213 - Icges(7,2) * t212) * t205;
t129 = -Icges(7,5) * t206 + (Icges(7,1) * t213 - Icges(7,4) * t212) * t205;
t56 = t127 * t296 + t128 * t159 + t129 * t160;
t5 = -t206 * t56 + (t227 * t28 + t230 * t29) * t205;
t30 = t161 * t91 + t162 * t93 + t295 * t89;
t31 = t161 * t92 + t162 * t94 + t295 * t90;
t57 = t127 * t295 + t128 * t161 + t129 * t162;
t6 = -t206 * t57 + (t227 * t30 + t230 * t31) * t205;
t321 = t6 * t295 + t5 * t296;
t320 = -t206 / 0.2e1;
t319 = t227 / 0.2e1;
t318 = -t230 / 0.2e1;
t317 = pkin(2) * t226;
t316 = pkin(4) * t206;
t224 = -qJ(3) - pkin(7);
t208 = t229 * pkin(2) + pkin(1);
t313 = rSges(3,1) * t229;
t311 = rSges(3,2) * t226;
t309 = t230 * rSges(3,3);
t40 = -t206 * t89 + (-t212 * t91 + t213 * t93) * t205;
t308 = t40 * t230;
t41 = -t206 * t90 + (-t212 * t92 + t213 * t94) * t205;
t307 = t41 * t227;
t49 = -t106 * t206 + (-t108 * t225 + t110 * t228) * t205;
t306 = t49 * t230;
t50 = -t107 * t206 + (-t109 * t225 + t111 * t228) * t205;
t305 = t50 * t227;
t117 = t205 * t213 * t129;
t293 = t212 * t128;
t63 = -t206 * t127 - t205 * t293 + t117;
t304 = t63 * t206;
t130 = -rSges(7,3) * t206 + (rSges(7,1) * t213 - rSges(7,2) * t212) * t205;
t70 = t130 * t296 + t206 * t95;
t302 = Icges(3,4) * t226;
t301 = Icges(3,4) * t229;
t300 = Icges(4,4) * t209;
t299 = Icges(4,4) * t210;
t132 = -Icges(6,6) * t206 + (Icges(6,4) * t228 - Icges(6,2) * t225) * t205;
t288 = t225 * t132;
t124 = t205 * t315 + t206 * t314;
t283 = -t124 - t130;
t134 = -rSges(6,3) * t206 + (rSges(6,1) * t228 - rSges(6,2) * t225) * t205;
t180 = t205 * pkin(4) - t206 * pkin(9);
t282 = -t134 - t180;
t241 = rSges(5,1) * t294 - rSges(5,2) * t295 + t227 * rSges(5,3);
t261 = rSges(5,1) * t206 - rSges(5,2) * t205;
t97 = t227 * (-rSges(5,3) * t230 + t227 * t261) + t230 * t241;
t201 = t230 * t208;
t218 = t230 * pkin(7);
t281 = t227 * (t218 + (-pkin(1) + t208) * t227) + t230 * (-pkin(1) * t230 + t201 - t327);
t280 = t221 * (pkin(9) * t205 + t316) + t230 * t279;
t278 = t227 * rSges(3,3) + t230 * t313;
t276 = t221 + t222;
t16 = t227 * t31 - t230 * t30;
t44 = t106 * t295 + t108 * t173 + t110 * t174;
t45 = t107 * t295 + t109 * t173 + t111 * t174;
t23 = t227 * t45 - t230 * t44;
t275 = (t221 * t142 + t16 + t23 + (t249 * t230 + (-t141 + t248) * t227) * t230) * t227;
t274 = -t180 + t283;
t273 = t296 / 0.2e1;
t272 = t295 / 0.2e1;
t271 = -rSges(4,1) * t209 - rSges(4,2) * t210 - t317;
t270 = (t40 + t56) * t273 + (t41 + t57) * t272;
t187 = pkin(3) * t210 + t208;
t181 = t230 * t187;
t219 = -pkin(8) + t224;
t269 = -t227 * t219 + t181;
t9 = -t304 + (t227 * t40 + t230 * t41) * t205;
t268 = -t206 * t9 + t321;
t267 = t230 * (t181 - t201) + t281 + (t187 - t208) * t221;
t266 = t15 * t273 + t16 * t272 + t5 * t318 + t6 * t319 + (t307 - t308) * t320;
t264 = -pkin(3) * t209 - t317;
t263 = -t311 + t313;
t258 = Icges(3,1) * t229 - t302;
t257 = Icges(4,1) * t210 - t300;
t255 = -Icges(3,2) * t226 + t301;
t254 = -Icges(4,2) * t209 + t299;
t177 = Icges(5,2) * t206 + t298;
t178 = Icges(5,1) * t205 + t297;
t243 = -t177 * t205 + t178 * t206;
t242 = t227 * rSges(4,3) + t230 * t262;
t179 = rSges(5,1) * t205 + rSges(5,2) * t206;
t240 = -t179 + t264;
t239 = -t180 + t264;
t237 = t267 + t280;
t236 = -t134 + t239;
t235 = t323 * t230 + t275;
t234 = t239 + t283;
t131 = -Icges(6,3) * t206 + (Icges(6,5) * t228 - Icges(6,6) * t225) * t205;
t133 = -Icges(6,5) * t206 + (Icges(6,1) * t228 - Icges(6,4) * t225) * t205;
t60 = t131 * t296 + t132 * t171 + t133 * t172;
t10 = -t206 * t60 + (t227 * t42 + t230 * t43) * t205;
t61 = t131 * t295 + t132 * t173 + t133 * t174;
t11 = -t206 * t61 + (t227 * t44 + t230 * t45) * t205;
t233 = t10 * t318 + t11 * t319 + t22 * t273 + t23 * t272 + (t305 - t306) * t320 + t266;
t176 = Icges(5,5) * t205 + Icges(5,6) * t206;
t232 = -t308 / 0.2e1 + t307 / 0.2e1 - t306 / 0.2e1 + t305 / 0.2e1 + (t144 * t206 + t146 * t205 + t176 * t227 + t230 * t243 + t57 + t61) * t319 + (t143 * t206 + t145 * t205 - t176 * t230 + t227 * t243 + t56 + t60) * t318;
t199 = rSges(2,1) * t230 - rSges(2,2) * t227;
t198 = -rSges(2,1) * t227 - rSges(2,2) * t230;
t197 = rSges(3,1) * t226 + rSges(3,2) * t229;
t148 = t271 * t230;
t147 = t271 * t227;
t140 = t327 + (pkin(1) - t311) * t230 + t278;
t139 = t309 + t218 + (-pkin(1) - t263) * t227;
t126 = t240 * t230;
t125 = t240 * t227;
t123 = -t224 * t227 + t201 + t242;
t122 = (rSges(4,3) - t224) * t230 + (-t208 - t262) * t227;
t121 = t205 * t228 * t133;
t116 = t230 * (-t230 * t311 + t278) + (t227 * t263 - t309) * t227;
t115 = t241 + t269;
t114 = (rSges(5,3) - t219) * t230 + (-t187 - t261) * t227;
t105 = t282 * t230;
t104 = t282 * t227;
t85 = t236 * t230;
t84 = t236 * t227;
t83 = t95 * t295;
t78 = t269 + t113 + t279;
t77 = -t230 * t219 + (-t316 - t187 + (-rSges(6,3) - pkin(9)) * t205) * t227 + t260;
t76 = t274 * t230;
t75 = t274 * t227;
t74 = -t113 * t206 - t134 * t295;
t73 = t112 * t206 + t134 * t296;
t72 = t230 * t242 + (-t230 * rSges(4,3) + t227 * t262) * t227 + t281;
t71 = -t130 * t295 - t206 * t96;
t69 = t234 * t230;
t68 = t234 * t227;
t67 = t269 + t331;
t66 = (pkin(5) * t225 - t219) * t230 + (-t206 * t207 - t187 + (-rSges(7,3) + t231) * t205) * t227 + t259;
t65 = -t206 * t131 - t205 * t288 + t121;
t64 = (t112 * t230 - t113 * t227) * t205;
t62 = -t296 * t96 + t83;
t55 = t280 + t326;
t46 = t267 + t97;
t35 = -t206 * t330 + t283 * t295;
t34 = t102 * t206 + t124 * t296 + t70;
t27 = t83 + (t102 * t230 - t227 * t330) * t205;
t26 = t280 + t322;
t25 = t237 + t326;
t18 = t237 + t322;
t1 = [t210 * (Icges(4,2) * t210 + t300) + t209 * (Icges(4,1) * t209 + t299) + t229 * (Icges(3,2) * t229 + t302) + t226 * (Icges(3,1) * t226 + t301) + Icges(2,3) + t117 + t121 + (-t127 - t131 + t177) * t206 + (t178 - t288 - t293) * t205 + m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t77 ^ 2 + t78 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2) + m(4) * (t122 ^ 2 + t123 ^ 2) + m(3) * (t139 ^ 2 + t140 ^ 2) + m(2) * (t198 ^ 2 + t199 ^ 2); t232 + m(3) * (-t139 * t230 - t140 * t227) * t197 + m(7) * (t66 * t69 + t67 * t68) + m(6) * (t77 * t85 + t78 * t84) + m(5) * (t114 * t126 + t115 * t125) + m(4) * (t122 * t148 + t123 * t147) + ((Icges(4,6) * t227 + t230 * t254) * t210 + (Icges(4,5) * t227 + t230 * t257) * t209 + (Icges(3,6) * t227 + t230 * t255) * t229 + (Icges(3,5) * t227 + t230 * t258) * t226) * t319 + ((-Icges(4,6) * t230 + t227 * t254) * t210 + (-Icges(4,5) * t230 + t227 * t257) * t209 + (-Icges(3,6) * t230 + t227 * t255) * t229 + (-Icges(3,5) * t230 + t227 * t258) * t226) * t318 + (Icges(3,5) * t226 + Icges(4,5) * t209 + Icges(3,6) * t229 + Icges(4,6) * t210) * (t222 / 0.2e1 + t221 / 0.2e1); m(7) * (t18 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t25 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t125 ^ 2 + t126 ^ 2 + t46 ^ 2) + m(4) * (t147 ^ 2 + t148 ^ 2 + t72 ^ 2) + m(3) * (t197 ^ 2 * t276 + t116 ^ 2) + t275 + t325 * t227 * t221 + (t324 * t222 + (t324 * t227 + t325 * t230) * t227 + t323) * t230; m(7) * (t227 * t66 - t230 * t67) + m(6) * (t227 * t77 - t230 * t78) + m(5) * (t114 * t227 - t115 * t230) + m(4) * (t122 * t227 - t123 * t230); m(7) * (t227 * t69 - t230 * t68) + m(6) * (t227 * t85 - t230 * t84) + m(5) * (-t125 * t230 + t126 * t227) + m(4) * (-t147 * t230 + t148 * t227); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t276; t232 + m(7) * (t66 * t76 + t67 * t75) + m(6) * (t104 * t78 + t105 * t77) + m(5) * (-t114 * t230 - t115 * t227) * t179; m(7) * (t26 * t18 + t68 * t75 + t69 * t76) + m(6) * (t104 * t84 + t105 * t85 + t55 * t25) + m(5) * (t46 * t97 + (-t125 * t227 - t126 * t230) * t179) + t235; m(6) * (-t104 * t230 + t105 * t227) + m(7) * (t227 * t76 - t230 * t75); m(7) * (t26 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t104 ^ 2 + t105 ^ 2 + t55 ^ 2) + m(5) * (t179 ^ 2 * t276 + t97 ^ 2) + t235; (-t63 - t65) * t206 + m(7) * (t34 * t66 + t35 * t67) + m(6) * (t73 * t77 + t74 * t78) + ((t61 / 0.2e1 + t50 / 0.2e1) * t230 + (t60 / 0.2e1 + t49 / 0.2e1) * t227) * t205 + t270; m(7) * (t18 * t27 + t34 * t69 + t35 * t68) + m(6) * (t64 * t25 + t73 * t85 + t74 * t84) + t233; m(6) * (t227 * t73 - t230 * t74) + m(7) * (t227 * t34 - t230 * t35); m(7) * (t27 * t26 + t34 * t76 + t35 * t75) + m(6) * (t104 * t74 + t105 * t73 + t64 * t55) + t233; (t65 * t206 - t9) * t206 + m(7) * (t27 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t64 ^ 2 + t73 ^ 2 + t74 ^ 2) + (t230 * t11 + t227 * t10 - t206 * (t227 * t49 + t230 * t50)) * t205 + t321; -t304 + m(7) * (t66 * t70 + t67 * t71) + t270; m(7) * (t62 * t18 + t68 * t71 + t69 * t70) + t266; m(7) * (t227 * t70 - t230 * t71); m(7) * (t62 * t26 + t70 * t76 + t71 * t75) + t266; m(7) * (t62 * t27 + t34 * t70 + t35 * t71) + t268; m(7) * (t62 ^ 2 + t70 ^ 2 + t71 ^ 2) + t268;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
