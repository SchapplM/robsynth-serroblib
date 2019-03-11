% Calculate joint inertia matrix for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:13
% EndTime: 2019-03-09 05:39:27
% DurationCPUTime: 5.85s
% Computational Cost: add. (36637->568), mult. (95497->789), div. (0->0), fcn. (125897->16), ass. (0->265)
t254 = sin(qJ(1));
t250 = cos(pkin(6));
t305 = cos(pkin(12));
t275 = t254 * t305;
t247 = sin(pkin(12));
t255 = cos(qJ(1));
t296 = t255 * t247;
t266 = t250 * t275 + t296;
t248 = sin(pkin(6));
t306 = cos(pkin(7));
t277 = t248 * t306;
t304 = sin(pkin(7));
t221 = t254 * t277 + t266 * t304;
t274 = t255 * t305;
t297 = t254 * t247;
t267 = -t250 * t274 + t297;
t220 = -t255 * t277 + t267 * t304;
t318 = m(7) / 0.2e1;
t319 = m(6) / 0.2e1;
t287 = t319 + t318;
t324 = 0.2e1 * t287;
t230 = -t250 * t297 + t274;
t253 = sin(qJ(3));
t260 = t266 * t306;
t276 = t248 * t304;
t311 = cos(qJ(3));
t214 = t230 * t311 + (t254 * t276 - t260) * t253;
t252 = sin(qJ(4));
t310 = cos(qJ(4));
t194 = t214 * t252 - t221 * t310;
t195 = t214 * t310 + t221 * t252;
t249 = cos(pkin(13));
t240 = pkin(5) * t249 + pkin(4);
t251 = -pkin(11) - qJ(5);
t271 = t311 * t304;
t268 = t248 * t271;
t213 = t230 * t253 - t254 * t268 + t260 * t311;
t246 = sin(pkin(13));
t302 = t213 * t246;
t245 = pkin(13) + qJ(6);
t241 = sin(t245);
t242 = cos(t245);
t150 = -t195 * t241 + t213 * t242;
t151 = t195 * t242 + t213 * t241;
t94 = t151 * rSges(7,1) + t150 * rSges(7,2) + t194 * rSges(7,3);
t323 = pkin(5) * t302 - t194 * t251 + t195 * t240 + t94;
t322 = m(3) / 0.2e1;
t321 = m(4) / 0.2e1;
t320 = m(5) / 0.2e1;
t229 = t250 * t296 + t275;
t262 = t267 * t306;
t212 = t229 * t311 + (-t255 * t276 - t262) * t253;
t192 = t212 * t252 - t220 * t310;
t270 = t306 * t305;
t219 = t250 * t304 * t253 + (t247 * t311 + t253 * t270) * t248;
t228 = t250 * t306 - t276 * t305;
t209 = t219 * t252 - t228 * t310;
t193 = t212 * t310 + t220 * t252;
t211 = t229 * t253 + t255 * t268 + t262 * t311;
t148 = -t193 * t241 + t211 * t242;
t149 = t193 * t242 + t211 * t241;
t87 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t192;
t89 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t192;
t91 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t192;
t28 = t148 * t89 + t149 * t91 + t192 * t87;
t88 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t194;
t90 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t194;
t92 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t194;
t29 = t148 * t90 + t149 * t92 + t192 * t88;
t210 = t219 * t310 + t228 * t252;
t300 = t247 * t248;
t218 = -t248 * t270 * t311 - t250 * t271 + t253 * t300;
t179 = -t210 * t241 + t218 * t242;
t180 = t210 * t242 + t218 * t241;
t116 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t209;
t117 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t209;
t118 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t209;
t44 = t116 * t192 + t117 * t148 + t118 * t149;
t1 = t192 * t28 + t194 * t29 + t209 * t44;
t317 = t1 / 0.2e1;
t30 = t150 * t89 + t151 * t91 + t194 * t87;
t31 = t150 * t90 + t151 * t92 + t194 * t88;
t45 = t116 * t194 + t117 * t150 + t118 * t151;
t2 = t192 * t30 + t194 * t31 + t209 * t45;
t316 = t2 / 0.2e1;
t37 = t179 * t89 + t180 * t91 + t209 * t87;
t38 = t179 * t90 + t180 * t92 + t209 * t88;
t54 = t209 * t116 + t179 * t117 + t180 * t118;
t50 = t54 * t209;
t11 = t37 * t192 + t38 * t194 + t50;
t315 = t11 / 0.2e1;
t314 = t192 / 0.2e1;
t313 = t194 / 0.2e1;
t312 = t209 / 0.2e1;
t309 = -pkin(4) + t240;
t186 = t192 * qJ(5);
t303 = t211 * t246;
t286 = pkin(5) * t303;
t269 = -t149 * rSges(7,1) - t148 * rSges(7,2);
t93 = rSges(7,3) * t192 - t269;
t308 = -t192 * t251 + t193 * t309 - t186 + t286 + t93;
t143 = t195 * pkin(4) + t194 * qJ(5);
t307 = -t143 + t323;
t301 = t218 * t246;
t299 = t248 * t254;
t298 = t248 * t255;
t154 = -t193 * t246 + t211 * t249;
t155 = t193 * t249 + t303;
t101 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t192;
t142 = pkin(4) * t193 + t186;
t295 = -t101 - t142;
t156 = -t195 * t246 + t213 * t249;
t157 = t195 * t249 + t302;
t102 = t157 * rSges(6,1) + t156 * rSges(6,2) + t194 * rSges(6,3);
t294 = -t102 - t143;
t119 = rSges(7,1) * t180 + rSges(7,2) * t179 + rSges(7,3) * t209;
t293 = pkin(5) * t301 + t309 * t210 + (-qJ(5) - t251) * t209 + t119;
t184 = -t210 * t246 + t218 * t249;
t185 = t210 * t249 + t301;
t124 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t209;
t175 = t210 * pkin(4) + t209 * qJ(5);
t292 = -t124 - t175;
t176 = t212 * pkin(3) + t211 * pkin(10);
t171 = t221 * t176;
t291 = t221 * t142 + t171;
t177 = t214 * pkin(3) + t213 * pkin(10);
t172 = t228 * t177;
t290 = t228 * t143 + t172;
t200 = pkin(3) * t219 + pkin(10) * t218;
t182 = t220 * t200;
t289 = t220 * t175 + t182;
t288 = t255 * pkin(1) + qJ(2) * t299;
t285 = -t142 - t308;
t284 = -t143 - t307;
t283 = t37 / 0.2e1 + t44 / 0.2e1;
t282 = t45 / 0.2e1 + t38 / 0.2e1;
t121 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t209;
t122 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t209;
t123 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t209;
t57 = t209 * t121 + t184 * t122 + t185 * t123;
t281 = -t175 - t293;
t158 = Icges(5,5) * t210 - Icges(5,6) * t209 + Icges(5,3) * t218;
t159 = Icges(5,4) * t210 - Icges(5,2) * t209 + Icges(5,6) * t218;
t160 = Icges(5,1) * t210 - Icges(5,4) * t209 + Icges(5,5) * t218;
t78 = t218 * t158 - t209 * t159 + t210 * t160;
t196 = Icges(4,5) * t219 - Icges(4,6) * t218 + Icges(4,3) * t228;
t197 = Icges(4,4) * t219 - Icges(4,2) * t218 + Icges(4,6) * t228;
t198 = Icges(4,1) * t219 - Icges(4,4) * t218 + Icges(4,5) * t228;
t280 = t228 * t196 - t218 * t197 + t219 * t198;
t132 = t195 * rSges(5,1) - t194 * rSges(5,2) + t213 * rSges(5,3);
t170 = t214 * rSges(4,1) - t213 * rSges(4,2) + t221 * rSges(4,3);
t278 = -t254 * pkin(1) + qJ(2) * t298;
t169 = rSges(4,1) * t212 - rSges(4,2) * t211 + rSges(4,3) * t220;
t131 = rSges(5,1) * t193 - rSges(5,2) * t192 + rSges(5,3) * t211;
t265 = -t229 * pkin(2) - pkin(9) * t220 + t278;
t95 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t192;
t97 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t192;
t99 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t192;
t39 = t184 * t97 + t185 * t99 + t209 * t95;
t46 = t121 * t192 + t122 * t154 + t123 * t155;
t125 = Icges(5,5) * t193 - Icges(5,6) * t192 + Icges(5,3) * t211;
t127 = Icges(5,4) * t193 - Icges(5,2) * t192 + Icges(5,6) * t211;
t129 = Icges(5,1) * t193 - Icges(5,4) * t192 + Icges(5,5) * t211;
t62 = t125 * t218 - t127 * t209 + t129 * t210;
t71 = t158 * t211 - t159 * t192 + t160 * t193;
t264 = t46 / 0.2e1 + t39 / 0.2e1 + t62 / 0.2e1 + t71 / 0.2e1 + t283;
t100 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t194;
t96 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t194;
t98 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t194;
t40 = t100 * t185 + t184 * t98 + t209 * t96;
t47 = t121 * t194 + t122 * t156 + t123 * t157;
t126 = Icges(5,5) * t195 - Icges(5,6) * t194 + Icges(5,3) * t213;
t128 = Icges(5,4) * t195 - Icges(5,2) * t194 + Icges(5,6) * t213;
t130 = Icges(5,1) * t195 - Icges(5,4) * t194 + Icges(5,5) * t213;
t63 = t126 * t218 - t128 * t209 + t130 * t210;
t72 = t158 * t213 - t159 * t194 + t160 * t195;
t263 = t63 / 0.2e1 + t40 / 0.2e1 + t72 / 0.2e1 + t47 / 0.2e1 + t282;
t258 = -t176 + t265;
t257 = t230 * pkin(2) + pkin(9) * t221 + t288;
t256 = t177 + t257;
t236 = rSges(2,1) * t255 - t254 * rSges(2,2);
t235 = -t254 * rSges(2,1) - rSges(2,2) * t255;
t208 = t230 * rSges(3,1) - rSges(3,2) * t266 + rSges(3,3) * t299 + t288;
t207 = -t229 * rSges(3,1) + rSges(3,2) * t267 + rSges(3,3) * t298 + t278;
t199 = rSges(4,1) * t219 - rSges(4,2) * t218 + rSges(4,3) * t228;
t168 = Icges(4,1) * t214 - Icges(4,4) * t213 + Icges(4,5) * t221;
t167 = Icges(4,1) * t212 - Icges(4,4) * t211 + Icges(4,5) * t220;
t166 = Icges(4,4) * t214 - Icges(4,2) * t213 + Icges(4,6) * t221;
t165 = Icges(4,4) * t212 - Icges(4,2) * t211 + Icges(4,6) * t220;
t164 = Icges(4,5) * t214 - Icges(4,6) * t213 + Icges(4,3) * t221;
t163 = Icges(4,5) * t212 - Icges(4,6) * t211 + Icges(4,3) * t220;
t161 = rSges(5,1) * t210 - rSges(5,2) * t209 + rSges(5,3) * t218;
t145 = t211 * t175;
t140 = t257 + t170;
t139 = -t169 + t265;
t135 = t218 * t143;
t133 = t213 * t142;
t114 = t170 * t228 - t199 * t221;
t113 = -t169 * t228 + t199 * t220;
t108 = t169 * t221 - t170 * t220;
t107 = t280 * t228;
t104 = t196 * t221 - t197 * t213 + t198 * t214;
t103 = t196 * t220 - t197 * t211 + t198 * t212;
t86 = t256 + t132;
t85 = -t131 + t258;
t82 = t132 * t218 - t161 * t213;
t81 = -t131 * t218 + t161 * t211;
t80 = t164 * t228 - t166 * t218 + t168 * t219;
t79 = t163 * t228 - t165 * t218 + t167 * t219;
t77 = t131 * t213 - t132 * t211;
t76 = t78 * t228;
t75 = t78 * t218;
t74 = t228 * t132 + t172 + (-t161 - t200) * t221;
t73 = t220 * t161 + t182 + (-t131 - t176) * t228;
t70 = t256 - t294;
t69 = t258 + t295;
t68 = -t119 * t194 + t209 * t94;
t67 = t119 * t192 - t209 * t93;
t66 = t221 * t131 + t171 + (-t132 - t177) * t220;
t65 = t256 + t323;
t64 = -t286 - t193 * t240 + (-rSges(7,3) + t251) * t192 + t258 + t269;
t61 = t126 * t213 - t128 * t194 + t130 * t195;
t60 = t125 * t213 - t127 * t194 + t129 * t195;
t59 = t126 * t211 - t128 * t192 + t130 * t193;
t58 = t125 * t211 - t127 * t192 + t129 * t193;
t56 = -t192 * t94 + t194 * t93;
t55 = t57 * t228;
t53 = t57 * t218;
t52 = t54 * t228;
t51 = t54 * t218;
t49 = t218 * t102 + t213 * t292 + t135;
t48 = t211 * t124 + t218 * t295 + t145;
t43 = t228 * t102 + (-t200 + t292) * t221 + t290;
t42 = t220 * t124 + (-t176 + t295) * t228 + t289;
t41 = t213 * t101 + t211 * t294 + t133;
t36 = t221 * t101 + (-t177 + t294) * t220 + t291;
t35 = t100 * t157 + t156 * t98 + t194 * t96;
t34 = t156 * t97 + t157 * t99 + t194 * t95;
t33 = t100 * t155 + t154 * t98 + t192 * t96;
t32 = t154 * t97 + t155 * t99 + t192 * t95;
t27 = t213 * t281 + t218 * t307 + t135;
t26 = t211 * t293 + t218 * t285 + t145;
t25 = t307 * t228 + (-t200 + t281) * t221 + t290;
t24 = t293 * t220 + (-t176 + t285) * t228 + t289;
t23 = t211 * t284 + t213 * t308 + t133;
t22 = t62 * t220 + t63 * t221 + t76;
t21 = t62 * t211 + t63 * t213 + t75;
t20 = t308 * t221 + (-t177 + t284) * t220 + t291;
t19 = t220 * t60 + t221 * t61 + t228 * t72;
t18 = t220 * t58 + t221 * t59 + t228 * t71;
t17 = t211 * t60 + t213 * t61 + t218 * t72;
t16 = t211 * t58 + t213 * t59 + t218 * t71;
t15 = t39 * t220 + t40 * t221 + t55;
t14 = t39 * t211 + t40 * t213 + t53;
t13 = t37 * t220 + t38 * t221 + t52;
t12 = t37 * t211 + t38 * t213 + t51;
t10 = t220 * t34 + t221 * t35 + t228 * t47;
t9 = t220 * t32 + t221 * t33 + t228 * t46;
t8 = t211 * t34 + t213 * t35 + t218 * t47;
t7 = t211 * t32 + t213 * t33 + t218 * t46;
t6 = t220 * t30 + t221 * t31 + t228 * t45;
t5 = t220 * t28 + t221 * t29 + t228 * t44;
t4 = t211 * t30 + t213 * t31 + t218 * t45;
t3 = t211 * t28 + t213 * t29 + t218 * t44;
t83 = [t248 * t305 * (Icges(3,6) * t250 + (Icges(3,4) * t247 + Icges(3,2) * t305) * t248) + t280 + t250 * (Icges(3,3) * t250 + (Icges(3,5) * t247 + Icges(3,6) * t305) * t248) + m(7) * (t64 ^ 2 + t65 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(3) * (t207 ^ 2 + t208 ^ 2) + m(2) * (t235 ^ 2 + t236 ^ 2) + Icges(2,3) + (Icges(3,5) * t250 + (Icges(3,1) * t247 + Icges(3,4) * t305) * t248) * t300 + t78 + t57 + t54; 0.2e1 * ((t254 * t64 - t255 * t65) * t318 + (t254 * t69 - t255 * t70) * t319 + (t254 * t85 - t255 * t86) * t320 + (t139 * t254 - t140 * t255) * t321 + (t207 * t254 - t208 * t255) * t322) * t248; 0.2e1 * (t322 + t321 + t320 + t287) * (t250 ^ 2 + (t254 ^ 2 + t255 ^ 2) * t248 ^ 2); t52 + t76 + t55 + t107 + m(7) * (t24 * t64 + t25 * t65) + m(6) * (t42 * t69 + t43 * t70) + m(5) * (t73 * t85 + t74 * t86) + m(4) * (t113 * t139 + t114 * t140) + (t80 / 0.2e1 + t104 / 0.2e1 + t263) * t221 + (t103 / 0.2e1 + t79 / 0.2e1 + t264) * t220; m(4) * (t108 * t250 + (t113 * t254 - t114 * t255) * t248) + m(5) * (t66 * t250 + (t254 * t73 - t255 * t74) * t248) + m(6) * (t36 * t250 + (t254 * t42 - t255 * t43) * t248) + m(7) * (t20 * t250 + (t24 * t254 - t25 * t255) * t248); (t107 + t13 + t15 + t22) * t228 + (t10 + t19 + t6 + (t164 * t221 - t166 * t213 + t168 * t214) * t221 + (t104 + t80) * t228) * t221 + (t5 + t18 + t9 + (t163 * t220 - t165 * t211 + t167 * t212) * t220 + (t103 + t79) * t228 + (t163 * t221 + t164 * t220 - t165 * t213 - t166 * t211 + t167 * t214 + t168 * t212) * t221) * t220 + m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t66 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t108 ^ 2 + t113 ^ 2 + t114 ^ 2); t51 + t75 + t53 + m(7) * (t26 * t64 + t27 * t65) + m(6) * (t48 * t69 + t49 * t70) + m(5) * (t81 * t85 + t82 * t86) + t263 * t213 + t264 * t211; m(5) * (t77 * t250 + (t254 * t81 - t255 * t82) * t248) + m(6) * (t41 * t250 + (t254 * t48 - t255 * t49) * t248) + m(7) * (t23 * t250 + (t254 * t26 - t255 * t27) * t248); (t12 / 0.2e1 + t14 / 0.2e1 + t21 / 0.2e1) * t228 + (t4 / 0.2e1 + t8 / 0.2e1 + t17 / 0.2e1) * t221 + (t3 / 0.2e1 + t7 / 0.2e1 + t16 / 0.2e1) * t220 + (t13 / 0.2e1 + t22 / 0.2e1 + t15 / 0.2e1) * t218 + (t6 / 0.2e1 + t19 / 0.2e1 + t10 / 0.2e1) * t213 + (t5 / 0.2e1 + t18 / 0.2e1 + t9 / 0.2e1) * t211 + m(7) * (t20 * t23 + t24 * t26 + t25 * t27) + m(6) * (t36 * t41 + t42 * t48 + t43 * t49) + m(5) * (t66 * t77 + t73 * t81 + t74 * t82); (t12 + t14 + t21) * t218 + (t4 + t17 + t8) * t213 + (t3 + t7 + t16) * t211 + m(7) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t77 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * (t192 * t65 + t194 * t64) + m(6) * (t192 * t70 + t194 * t69); (t209 * t250 + (-t192 * t255 + t194 * t254) * t248) * t324; m(7) * (t192 * t25 + t194 * t24 + t20 * t209) + m(6) * (t192 * t43 + t194 * t42 + t209 * t36); m(7) * (t192 * t27 + t194 * t26 + t209 * t23) + m(6) * (t192 * t49 + t194 * t48 + t209 * t41); (t192 ^ 2 + t194 ^ 2 + t209 ^ 2) * t324; t50 + m(7) * (t64 * t67 + t65 * t68) + t282 * t194 + t283 * t192; m(7) * (t56 * t250 + (t254 * t67 - t255 * t68) * t248); m(7) * (t20 * t56 + t24 * t67 + t25 * t68) + t221 * t316 + t220 * t317 + t5 * t314 + t6 * t313 + t13 * t312 + t228 * t315; m(7) * (t23 * t56 + t26 * t67 + t27 * t68) + t4 * t313 + t3 * t314 + t211 * t317 + t213 * t316 + t12 * t312 + t218 * t315; m(7) * (t192 * t68 + t194 * t67 + t209 * t56); t194 * t2 + t192 * t1 + t209 * t11 + m(7) * (t56 ^ 2 + t67 ^ 2 + t68 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t83(1) t83(2) t83(4) t83(7) t83(11) t83(16); t83(2) t83(3) t83(5) t83(8) t83(12) t83(17); t83(4) t83(5) t83(6) t83(9) t83(13) t83(18); t83(7) t83(8) t83(9) t83(10) t83(14) t83(19); t83(11) t83(12) t83(13) t83(14) t83(15) t83(20); t83(16) t83(17) t83(18) t83(19) t83(20) t83(21);];
Mq  = res;
