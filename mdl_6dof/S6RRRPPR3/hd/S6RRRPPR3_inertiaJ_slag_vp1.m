% Calculate joint inertia matrix for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:26
% EndTime: 2019-03-09 15:28:35
% DurationCPUTime: 4.83s
% Computational Cost: add. (5068->352), mult. (6411->523), div. (0->0), fcn. (6477->8), ass. (0->173)
t301 = Icges(4,4) + Icges(6,4) - Icges(5,5);
t300 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t299 = -Icges(6,1) - Icges(4,2) - Icges(5,3);
t184 = qJ(2) + qJ(3);
t176 = sin(t184);
t298 = t301 * t176;
t177 = cos(t184);
t297 = t301 * t177;
t296 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t295 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t294 = t176 * t299 + t297;
t293 = -t177 * t300 + t298;
t292 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t187 = sin(qJ(1));
t190 = cos(qJ(1));
t291 = t294 * t187 - t295 * t190;
t290 = t295 * t187 + t294 * t190;
t289 = t293 * t187 + t296 * t190;
t288 = -t296 * t187 + t293 * t190;
t287 = t295 * t176 - t296 * t177;
t286 = t177 * t299 - t298;
t285 = t176 * t300 + t297;
t284 = t292 * t187 - t287 * t190;
t283 = t287 * t187 + t292 * t190;
t282 = t291 * t176 + t289 * t177;
t281 = t290 * t176 + t288 * t177;
t182 = t187 ^ 2;
t280 = t187 * pkin(7);
t249 = t177 * t190;
t251 = t176 * t190;
t185 = sin(qJ(6));
t246 = t190 * t185;
t188 = cos(qJ(6));
t247 = t187 * t188;
t137 = -t176 * t246 - t247;
t245 = t190 * t188;
t248 = t187 * t185;
t138 = t176 * t245 - t248;
t74 = t138 * rSges(7,1) + t137 * rSges(7,2) + rSges(7,3) * t249;
t279 = pkin(5) * t251 + pkin(9) * t249 + t74;
t278 = -t296 * t176 - t295 * t177;
t277 = t176 * t286 + t285 * t177;
t99 = t176 * rSges(7,3) + (-rSges(7,1) * t188 + rSges(7,2) * t185) * t177;
t276 = t177 * pkin(5) - t176 * pkin(9) - t99;
t183 = t190 ^ 2;
t135 = -t176 * t248 + t245;
t136 = t176 * t247 + t246;
t250 = t177 * t187;
t67 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t250;
t69 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t250;
t71 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t250;
t16 = t135 * t69 + t136 * t71 + t67 * t250;
t68 = Icges(7,5) * t138 + Icges(7,6) * t137 + Icges(7,3) * t249;
t70 = Icges(7,4) * t138 + Icges(7,2) * t137 + Icges(7,6) * t249;
t72 = Icges(7,1) * t138 + Icges(7,4) * t137 + Icges(7,5) * t249;
t17 = t135 * t70 + t136 * t72 + t68 * t250;
t8 = -t16 * t190 + t17 * t187;
t275 = -t8 + t283 * t183 + (t281 * t187 + (-t282 + t284) * t190) * t187;
t274 = m(5) / 0.2e1;
t273 = m(6) / 0.2e1;
t272 = m(7) / 0.2e1;
t271 = -pkin(3) - pkin(4);
t270 = t187 / 0.2e1;
t269 = -t190 / 0.2e1;
t186 = sin(qJ(2));
t268 = pkin(2) * t186;
t96 = Icges(7,3) * t176 + (-Icges(7,5) * t188 + Icges(7,6) * t185) * t177;
t97 = Icges(7,6) * t176 + (-Icges(7,4) * t188 + Icges(7,2) * t185) * t177;
t267 = t177 * t185 * t97 + t176 * t96;
t189 = cos(qJ(2));
t174 = t189 * pkin(2) + pkin(1);
t162 = t190 * t174;
t181 = t190 * pkin(7);
t266 = t187 * (t181 + (-pkin(1) + t174) * t187) + t190 * (-t190 * pkin(1) + t162 - t280);
t197 = rSges(4,1) * t249 - rSges(4,2) * t251 + t187 * rSges(4,3);
t223 = rSges(4,1) * t177 - rSges(4,2) * t176;
t62 = t187 * (-t190 * rSges(4,3) + t223 * t187) + t190 * t197;
t265 = rSges(3,1) * t189;
t264 = rSges(3,2) * t186;
t98 = Icges(7,5) * t176 + (-Icges(7,1) * t188 + Icges(7,4) * t185) * t177;
t263 = t188 * t98;
t262 = t190 * rSges(3,3);
t26 = t176 * t67 + (t185 * t69 - t188 * t71) * t177;
t261 = t26 * t190;
t27 = t176 * t68 + (t185 * t70 - t188 * t72) * t177;
t260 = t27 * t187;
t259 = Icges(3,4) * t186;
t258 = Icges(3,4) * t189;
t191 = -pkin(8) - pkin(7);
t244 = -qJ(5) - t191;
t240 = pkin(3) * t249 + qJ(4) * t251;
t243 = t182 * (pkin(3) * t177 + qJ(4) * t176) + t190 * t240;
t149 = t176 * pkin(3) - t177 * qJ(4);
t242 = -t176 * rSges(5,1) + t177 * rSges(5,3) - t149;
t239 = t187 * rSges(3,3) + t190 * t265;
t238 = t182 + t183;
t237 = t273 + t272;
t236 = -rSges(6,3) + t244;
t235 = rSges(5,1) * t249 + t187 * rSges(5,2) + rSges(5,3) * t251;
t18 = t137 * t69 + t138 * t71 + t67 * t249;
t19 = t137 * t70 + t138 * t72 + t68 * t249;
t9 = -t18 * t190 + t19 * t187;
t234 = (t9 + t284 * t182 + ((-t281 + t283) * t187 + t282 * t190) * t190) * t187;
t151 = t176 * rSges(4,1) + t177 * rSges(4,2);
t233 = -t151 - t268;
t232 = -pkin(4) * t176 - t149;
t231 = -t187 * t191 + t162;
t39 = t187 * (-t190 * rSges(5,2) + (rSges(5,1) * t177 + rSges(5,3) * t176) * t187) + t190 * t235 + t243;
t32 = t135 * t97 + t136 * t98 + t96 * t250;
t3 = t32 * t176 + (t16 * t187 + t17 * t190) * t177;
t33 = t137 * t97 + t138 * t98 + t96 * t249;
t4 = t33 * t176 + (t18 * t187 + t19 * t190) * t177;
t230 = t3 * t269 + t4 * t270 + t176 * (t260 - t261) / 0.2e1 + t8 * t250 / 0.2e1 + t9 * t249 / 0.2e1;
t170 = pkin(4) * t249;
t229 = t187 * pkin(4) * t250 + t190 * t170 + t243;
t228 = t162 + t170 + t240;
t227 = t242 - t268;
t152 = -t177 * rSges(6,1) - t176 * rSges(6,2);
t226 = -t152 + t232;
t225 = rSges(6,1) * t251 - rSges(6,2) * t249;
t224 = -t264 + t265;
t222 = -t136 * rSges(7,1) - t135 * rSges(7,2);
t221 = Icges(3,1) * t189 - t259;
t217 = -Icges(3,2) * t186 + t258;
t213 = Icges(3,5) * t189 - Icges(3,6) * t186;
t198 = t232 + t276;
t196 = t232 - t268;
t29 = (rSges(6,1) * t176 - rSges(6,2) * t177) * t182 + t190 * t225 + t229;
t195 = -t152 + t196;
t194 = t196 + t276;
t73 = rSges(7,3) * t250 - t222;
t12 = t187 * t73 + t182 * (pkin(5) * t176 + pkin(9) * t177) + t229 + t279 * t190;
t193 = t275 * t190 + t234;
t192 = t260 / 0.2e1 - t261 / 0.2e1 + (-t288 * t176 + t290 * t177 - t278 * t187 + t277 * t190 + t33) * t270 + (-t289 * t176 + t291 * t177 + t277 * t187 + t278 * t190 + t32) * t269;
t160 = t190 * rSges(2,1) - t187 * rSges(2,2);
t159 = -t187 * rSges(2,1) - t190 * rSges(2,2);
t158 = t186 * rSges(3,1) + t189 * rSges(3,2);
t126 = Icges(3,3) * t187 + t213 * t190;
t125 = -Icges(3,3) * t190 + t213 * t187;
t101 = t233 * t190;
t100 = t233 * t187;
t87 = t280 + (pkin(1) - t264) * t190 + t239;
t86 = t262 + t181 + (-pkin(1) - t224) * t187;
t83 = t242 * t190;
t82 = t242 * t187;
t81 = t197 + t231;
t80 = (rSges(4,3) - t191) * t190 + (-t174 - t223) * t187;
t79 = t227 * t190;
t78 = t227 * t187;
t77 = t226 * t190;
t76 = t226 * t187;
t75 = t190 * (-t190 * t264 + t239) + (t224 * t187 - t262) * t187;
t66 = t195 * t190;
t65 = t195 * t187;
t61 = t231 + t235 + t240;
t60 = (rSges(5,2) - t191) * t190 + (-t174 + (-rSges(5,1) - pkin(3)) * t177 + (-rSges(5,3) - qJ(4)) * t176) * t187;
t47 = t236 * t187 + t225 + t228;
t46 = t236 * t190 + (-t174 + (-rSges(6,1) - qJ(4)) * t176 + (rSges(6,2) + t271) * t177) * t187;
t45 = t198 * t190;
t44 = t198 * t187;
t43 = t194 * t190;
t42 = t194 * t187;
t41 = t176 * t74 - t99 * t249;
t40 = -t176 * t73 + t99 * t250;
t38 = (-t177 * t263 + t267) * t176;
t37 = t62 + t266;
t36 = t244 * t187 + t228 + t279;
t35 = t244 * t190 + (-t174 + (-pkin(5) - qJ(4)) * t176 + (-rSges(7,3) - pkin(9) + t271) * t177) * t187 + t222;
t34 = (-t187 * t74 + t190 * t73) * t177;
t28 = t39 + t266;
t23 = t29 + t266;
t11 = t12 + t266;
t1 = [t189 * (Icges(3,2) * t189 + t259) + t186 * (Icges(3,1) * t186 + t258) + Icges(2,3) + t285 * t176 + (-t263 - t286) * t177 + m(7) * (t35 ^ 2 + t36 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t80 ^ 2 + t81 ^ 2) + m(3) * (t86 ^ 2 + t87 ^ 2) + m(2) * (t159 ^ 2 + t160 ^ 2) + t267; (t183 / 0.2e1 + t182 / 0.2e1) * (Icges(3,5) * t186 + Icges(3,6) * t189) + t192 + m(4) * (t100 * t81 + t101 * t80) + m(5) * (t79 * t60 + t78 * t61) + m(6) * (t66 * t46 + t65 * t47) + m(7) * (t35 * t43 + t42 * t36) + m(3) * (-t187 * t87 - t190 * t86) * t158 + (t189 * (Icges(3,6) * t187 + t217 * t190) + t186 * (Icges(3,5) * t187 + t221 * t190)) * t270 + (t189 * (-Icges(3,6) * t190 + t217 * t187) + t186 * (-Icges(3,5) * t190 + t221 * t187)) * t269; m(7) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t23 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t28 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t37 ^ 2) + m(3) * (t238 * t158 ^ 2 + t75 ^ 2) + t187 * t182 * t126 + t234 + (-t183 * t125 + (-t187 * t125 + t190 * t126) * t187 + t275) * t190; t192 + m(7) * (t35 * t45 + t36 * t44) + m(5) * (t83 * t60 + t82 * t61) + m(6) * (t77 * t46 + t76 * t47) + m(4) * (-t187 * t81 - t190 * t80) * t151; m(7) * (t12 * t11 + t42 * t44 + t43 * t45) + m(6) * (t29 * t23 + t76 * t65 + t77 * t66) + m(5) * (t39 * t28 + t82 * t78 + t83 * t79) + m(4) * (t62 * t37 + (-t100 * t187 - t101 * t190) * t151) + t193; m(7) * (t12 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t29 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t39 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(4) * (t238 * t151 ^ 2 + t62 ^ 2) + t193; 0.2e1 * ((t187 * t36 + t190 * t35) * t272 + (t187 * t61 + t190 * t60) * t274 + (t187 * t47 + t190 * t46) * t273) * t176; m(7) * (-t177 * t11 + (t187 * t42 + t190 * t43) * t176) + m(6) * (-t177 * t23 + (t187 * t65 + t190 * t66) * t176) + m(5) * (-t177 * t28 + (t187 * t78 + t190 * t79) * t176); m(7) * (-t177 * t12 + (t187 * t44 + t190 * t45) * t176) + m(6) * (-t177 * t29 + (t187 * t76 + t190 * t77) * t176) + m(5) * (-t177 * t39 + (t187 * t82 + t190 * t83) * t176); 0.2e1 * (t274 + t237) * (t238 * t176 ^ 2 + t177 ^ 2); m(7) * (-t187 * t35 + t190 * t36) + m(6) * (-t187 * t46 + t190 * t47); m(7) * (-t187 * t43 + t190 * t42) + m(6) * (-t187 * t66 + t190 * t65); m(7) * (-t187 * t45 + t190 * t44) + m(6) * (-t187 * t77 + t190 * t76); 0; 0.2e1 * t237 * t238; t38 + m(7) * (t35 * t40 + t36 * t41) + ((t27 / 0.2e1 + t33 / 0.2e1) * t190 + (t26 / 0.2e1 + t32 / 0.2e1) * t187) * t177; m(7) * (t34 * t11 + t40 * t43 + t41 * t42) + t230; m(7) * (t34 * t12 + t40 * t45 + t41 * t44) + t230; m(7) * (-t34 * t177 + (t187 * t41 + t190 * t40) * t176); m(7) * (-t40 * t187 + t41 * t190); t176 * t38 + m(7) * (t34 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t190 * t4 + t187 * t3 + t176 * (t187 * t26 + t190 * t27)) * t177;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
