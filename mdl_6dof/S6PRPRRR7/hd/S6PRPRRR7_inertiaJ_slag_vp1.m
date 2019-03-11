% Calculate joint inertia matrix for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:50
% EndTime: 2019-03-08 20:51:12
% DurationCPUTime: 9.85s
% Computational Cost: add. (91926->598), mult. (261548->879), div. (0->0), fcn. (348872->18), ass. (0->270)
t226 = sin(pkin(13));
t232 = cos(pkin(6));
t236 = sin(qJ(2));
t230 = cos(pkin(13));
t238 = cos(qJ(2));
t277 = t238 * t230;
t221 = -t226 * t236 + t232 * t277;
t227 = sin(pkin(7));
t228 = sin(pkin(6));
t231 = cos(pkin(7));
t282 = t228 * t231;
t211 = -t221 * t227 - t230 * t282;
t278 = t238 * t226;
t279 = t232 * t236;
t222 = t230 * t279 + t278;
t196 = pkin(2) * t222 + qJ(3) * t211;
t223 = -t230 * t236 - t232 * t278;
t212 = -t223 * t227 + t226 * t282;
t224 = -t226 * t279 + t277;
t197 = pkin(2) * t224 + qJ(3) * t212;
t283 = t228 * t230;
t286 = t226 * t228;
t271 = t196 * t286 + t197 * t283;
t225 = sin(pkin(14));
t229 = cos(pkin(14));
t285 = t227 * t228;
t244 = t223 * t231 + t226 * t285;
t193 = -t224 * t225 + t229 * t244;
t194 = t224 * t229 + t225 * t244;
t235 = sin(qJ(4));
t287 = sin(pkin(8));
t295 = cos(qJ(4));
t249 = t295 * t287;
t288 = cos(pkin(8));
t250 = t288 * t295;
t155 = -t193 * t250 + t194 * t235 - t212 * t249;
t156 = t194 * t295 + (t193 * t288 + t212 * t287) * t235;
t174 = -t193 * t287 + t212 * t288;
t120 = rSges(5,1) * t156 - rSges(5,2) * t155 + rSges(5,3) * t174;
t234 = sin(qJ(5));
t294 = cos(qJ(5));
t140 = t156 * t234 - t174 * t294;
t141 = t156 * t294 + t174 * t234;
t101 = rSges(6,1) * t141 - rSges(6,2) * t140 + rSges(6,3) * t155;
t233 = sin(qJ(6));
t237 = cos(qJ(6));
t106 = -t141 * t233 + t155 * t237;
t107 = t141 * t237 + t155 * t233;
t83 = rSges(7,1) * t107 + rSges(7,2) * t106 + rSges(7,3) * t140;
t290 = pkin(5) * t141 + pkin(12) * t140 + t83;
t240 = -m(6) * t101 - m(7) * t290;
t307 = -m(5) * t120 + t240;
t306 = m(6) / 0.2e1 + m(7) / 0.2e1;
t305 = m(4) / 0.2e1;
t304 = m(5) / 0.2e1;
t245 = t221 * t231 - t227 * t283;
t191 = -t222 * t225 + t229 * t245;
t192 = t222 * t229 + t225 * t245;
t154 = t192 * t295 + (t191 * t288 + t211 * t287) * t235;
t173 = -t191 * t287 + t211 * t288;
t138 = t154 * t234 - t173 * t294;
t280 = t231 * t238;
t284 = t227 * t232;
t209 = t229 * t284 + (-t225 * t236 + t229 * t280) * t228;
t281 = t228 * t236;
t210 = t229 * t281 + (t228 * t280 + t284) * t225;
t219 = t231 * t232 - t238 * t285;
t172 = t210 * t295 + (t209 * t288 + t219 * t287) * t235;
t190 = -t209 * t287 + t219 * t288;
t149 = t172 * t234 - t190 * t294;
t139 = t154 * t294 + t173 * t234;
t153 = -t191 * t250 + t192 * t235 - t211 * t249;
t104 = -t139 * t233 + t153 * t237;
t105 = t139 * t237 + t153 * t233;
t76 = Icges(7,5) * t105 + Icges(7,6) * t104 + Icges(7,3) * t138;
t78 = Icges(7,4) * t105 + Icges(7,2) * t104 + Icges(7,6) * t138;
t80 = Icges(7,1) * t105 + Icges(7,4) * t104 + Icges(7,5) * t138;
t27 = t104 * t78 + t105 * t80 + t138 * t76;
t77 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t140;
t79 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t140;
t81 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t140;
t28 = t104 * t79 + t105 * t81 + t138 * t77;
t150 = t172 * t294 + t190 * t234;
t171 = -t209 * t250 + t210 * t235 - t219 * t249;
t136 = -t150 * t233 + t171 * t237;
t137 = t150 * t237 + t171 * t233;
t90 = Icges(7,5) * t137 + Icges(7,6) * t136 + Icges(7,3) * t149;
t91 = Icges(7,4) * t137 + Icges(7,2) * t136 + Icges(7,6) * t149;
t92 = Icges(7,1) * t137 + Icges(7,4) * t136 + Icges(7,5) * t149;
t41 = t104 * t91 + t105 * t92 + t138 * t90;
t1 = t138 * t27 + t140 * t28 + t149 * t41;
t301 = t1 / 0.2e1;
t29 = t106 * t78 + t107 * t80 + t140 * t76;
t30 = t106 * t79 + t107 * t81 + t140 * t77;
t42 = t106 * t91 + t107 * t92 + t140 * t90;
t2 = t138 * t29 + t140 * t30 + t149 * t42;
t300 = t2 / 0.2e1;
t34 = t136 * t78 + t137 * t80 + t149 * t76;
t35 = t136 * t79 + t137 * t81 + t149 * t77;
t45 = t136 * t91 + t137 * t92 + t149 * t90;
t9 = t138 * t34 + t140 * t35 + t149 * t45;
t299 = t9 / 0.2e1;
t298 = t138 / 0.2e1;
t297 = t140 / 0.2e1;
t296 = t149 / 0.2e1;
t82 = rSges(7,1) * t105 + rSges(7,2) * t104 + rSges(7,3) * t138;
t291 = pkin(5) * t139 + pkin(12) * t138 + t82;
t93 = rSges(7,1) * t137 + rSges(7,2) * t136 + rSges(7,3) * t149;
t289 = pkin(5) * t150 + pkin(12) * t149 + t93;
t159 = pkin(3) * t194 + pkin(10) * t174;
t195 = t232 * t197;
t275 = t159 * t232 + t195;
t158 = pkin(3) * t192 + pkin(10) * t173;
t274 = -t158 - t196;
t213 = pkin(2) * t281 + qJ(3) * t219;
t273 = -pkin(3) * t210 - pkin(10) * t190 - t213;
t272 = 0.2e1 * t271;
t94 = Icges(6,5) * t139 - Icges(6,6) * t138 + Icges(6,3) * t153;
t96 = Icges(6,4) * t139 - Icges(6,2) * t138 + Icges(6,6) * t153;
t98 = Icges(6,1) * t139 - Icges(6,4) * t138 + Icges(6,5) * t153;
t47 = -t138 * t96 + t139 * t98 + t153 * t94;
t95 = Icges(6,5) * t141 - Icges(6,6) * t140 + Icges(6,3) * t155;
t97 = Icges(6,4) * t141 - Icges(6,2) * t140 + Icges(6,6) * t155;
t99 = Icges(6,1) * t141 - Icges(6,4) * t140 + Icges(6,5) * t155;
t48 = -t138 * t97 + t139 * t99 + t153 * t95;
t109 = Icges(6,5) * t150 - Icges(6,6) * t149 + Icges(6,3) * t171;
t110 = Icges(6,4) * t150 - Icges(6,2) * t149 + Icges(6,6) * t171;
t111 = Icges(6,1) * t150 - Icges(6,4) * t149 + Icges(6,5) * t171;
t59 = t109 * t153 - t110 * t138 + t111 * t139;
t13 = t153 * t47 + t155 * t48 + t171 * t59;
t3 = t153 * t27 + t155 * t28 + t171 * t41;
t269 = t3 / 0.2e1 + t13 / 0.2e1;
t49 = -t140 * t96 + t141 * t98 + t155 * t94;
t50 = -t140 * t97 + t141 * t99 + t155 * t95;
t60 = t109 * t155 - t110 * t140 + t111 * t141;
t14 = t153 * t49 + t155 * t50 + t171 * t60;
t4 = t153 * t29 + t155 * t30 + t171 * t42;
t268 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t173 * t47 + t174 * t48 + t190 * t59;
t5 = t173 * t27 + t174 * t28 + t190 * t41;
t267 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t173 * t49 + t174 * t50 + t190 * t60;
t6 = t173 * t29 + t174 * t30 + t190 * t42;
t266 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t232 * t59 + (t226 * t48 - t230 * t47) * t228;
t7 = t232 * t41 + (t226 * t28 - t230 * t27) * t228;
t265 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t232 * t60 + (t226 * t50 - t230 * t49) * t228;
t8 = t232 * t42 + (t226 * t30 - t230 * t29) * t228;
t264 = t8 / 0.2e1 + t18 / 0.2e1;
t10 = t153 * t34 + t155 * t35 + t171 * t45;
t51 = -t149 * t96 + t150 * t98 + t171 * t94;
t52 = -t149 * t97 + t150 * t99 + t171 * t95;
t63 = t109 * t171 - t110 * t149 + t111 * t150;
t19 = t153 * t51 + t155 * t52 + t171 * t63;
t263 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t173 * t34 + t174 * t35 + t190 * t45;
t20 = t173 * t51 + t174 * t52 + t190 * t63;
t262 = t11 / 0.2e1 + t20 / 0.2e1;
t12 = t232 * t45 + (t226 * t35 - t230 * t34) * t228;
t21 = t232 * t63 + (t226 * t52 - t230 * t51) * t228;
t261 = t12 / 0.2e1 + t21 / 0.2e1;
t134 = pkin(4) * t156 + pkin(11) * t155;
t260 = t134 * t232 + t275;
t133 = pkin(4) * t154 + pkin(11) * t153;
t259 = -t133 + t274;
t146 = pkin(4) * t172 + pkin(11) * t171;
t258 = -t146 + t273;
t254 = (-rSges(4,1) * t210 - rSges(4,2) * t209 - rSges(4,3) * t219 - t213) * t228;
t151 = t158 * t286;
t152 = t159 * t283;
t253 = 0.2e1 * t151 + 0.2e1 * t152 + t272;
t252 = t151 + t152 + t271;
t145 = rSges(5,1) * t172 - rSges(5,2) * t171 + rSges(5,3) * t190;
t251 = (-t145 + t273) * t228;
t46 = -t138 * t83 + t140 * t82;
t112 = rSges(6,1) * t150 - rSges(6,2) * t149 + rSges(6,3) * t171;
t248 = (-t112 + t258) * t228;
t129 = t133 * t286;
t130 = t134 * t283;
t246 = t129 + t130 + t252;
t243 = (t258 - t289) * t228;
t242 = 0.2e1 * t305 + 0.2e1 * t304 + 0.2e1 * t306;
t100 = rSges(6,1) * t139 - rSges(6,2) * t138 + rSges(6,3) * t153;
t241 = m(6) * t100 + m(7) * t291;
t119 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t173;
t239 = m(5) * t119 + t241;
t218 = t232 * rSges(3,3) + (rSges(3,1) * t236 + rSges(3,2) * t238) * t228;
t217 = Icges(3,5) * t232 + (Icges(3,1) * t236 + Icges(3,4) * t238) * t228;
t216 = Icges(3,6) * t232 + (Icges(3,4) * t236 + Icges(3,2) * t238) * t228;
t215 = Icges(3,3) * t232 + (Icges(3,5) * t236 + Icges(3,6) * t238) * t228;
t205 = rSges(3,1) * t224 + rSges(3,2) * t223 + rSges(3,3) * t286;
t204 = rSges(3,1) * t222 + rSges(3,2) * t221 - rSges(3,3) * t283;
t203 = Icges(3,1) * t224 + Icges(3,4) * t223 + Icges(3,5) * t286;
t202 = Icges(3,1) * t222 + Icges(3,4) * t221 - Icges(3,5) * t283;
t201 = Icges(3,4) * t224 + Icges(3,2) * t223 + Icges(3,6) * t286;
t200 = Icges(3,4) * t222 + Icges(3,2) * t221 - Icges(3,6) * t283;
t199 = Icges(3,5) * t224 + Icges(3,6) * t223 + Icges(3,3) * t286;
t198 = Icges(3,5) * t222 + Icges(3,6) * t221 - Icges(3,3) * t283;
t182 = -t204 * t232 - t218 * t283;
t181 = t205 * t232 - t218 * t286;
t178 = Icges(4,1) * t210 + Icges(4,4) * t209 + Icges(4,5) * t219;
t177 = Icges(4,4) * t210 + Icges(4,2) * t209 + Icges(4,6) * t219;
t176 = Icges(4,5) * t210 + Icges(4,6) * t209 + Icges(4,3) * t219;
t168 = (t204 * t226 + t205 * t230) * t228;
t167 = rSges(4,1) * t194 + rSges(4,2) * t193 + rSges(4,3) * t212;
t166 = rSges(4,1) * t192 + rSges(4,2) * t191 + rSges(4,3) * t211;
t165 = Icges(4,1) * t194 + Icges(4,4) * t193 + Icges(4,5) * t212;
t164 = Icges(4,1) * t192 + Icges(4,4) * t191 + Icges(4,5) * t211;
t163 = Icges(4,4) * t194 + Icges(4,2) * t193 + Icges(4,6) * t212;
t162 = Icges(4,4) * t192 + Icges(4,2) * t191 + Icges(4,6) * t211;
t161 = Icges(4,5) * t194 + Icges(4,6) * t193 + Icges(4,3) * t212;
t160 = Icges(4,5) * t192 + Icges(4,6) * t191 + Icges(4,3) * t211;
t144 = Icges(5,1) * t172 - Icges(5,4) * t171 + Icges(5,5) * t190;
t143 = Icges(5,4) * t172 - Icges(5,2) * t171 + Icges(5,6) * t190;
t142 = Icges(5,5) * t172 - Icges(5,6) * t171 + Icges(5,3) * t190;
t135 = t173 * t146;
t126 = (-t166 - t196) * t232 + t230 * t254;
t125 = t167 * t232 + t226 * t254 + t195;
t124 = t190 * t134;
t123 = t174 * t133;
t118 = Icges(5,1) * t156 - Icges(5,4) * t155 + Icges(5,5) * t174;
t117 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t173;
t116 = Icges(5,4) * t156 - Icges(5,2) * t155 + Icges(5,6) * t174;
t115 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t173;
t114 = Icges(5,5) * t156 - Icges(5,6) * t155 + Icges(5,3) * t174;
t113 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t173;
t108 = (t166 * t226 + t167 * t230) * t228 + t271;
t89 = t120 * t190 - t145 * t174;
t88 = -t119 * t190 + t145 * t173;
t87 = t142 * t190 - t143 * t171 + t144 * t172;
t86 = t119 * t174 - t120 * t173;
t85 = (-t119 + t274) * t232 + t230 * t251;
t84 = t120 * t232 + t226 * t251 + t275;
t75 = t142 * t174 - t143 * t155 + t144 * t156;
t74 = t142 * t173 - t143 * t153 + t144 * t154;
t73 = (t119 * t226 + t120 * t230) * t228 + t252;
t72 = t101 * t171 - t112 * t155;
t71 = -t100 * t171 + t112 * t153;
t70 = t114 * t190 - t116 * t171 + t118 * t172;
t69 = t113 * t190 - t115 * t171 + t117 * t172;
t68 = t114 * t174 - t116 * t155 + t118 * t156;
t67 = t113 * t174 - t115 * t155 + t117 * t156;
t66 = t114 * t173 - t116 * t153 + t118 * t154;
t65 = t113 * t173 - t115 * t153 + t117 * t154;
t64 = t100 * t155 - t101 * t153;
t62 = t101 * t190 + t124 + (-t112 - t146) * t174;
t61 = t112 * t173 + t135 + (-t100 - t133) * t190;
t58 = (-t100 + t259) * t232 + t230 * t248;
t57 = t101 * t232 + t226 * t248 + t260;
t56 = -t140 * t93 + t149 * t83;
t55 = t138 * t93 - t149 * t82;
t54 = t100 * t174 + t123 + (-t101 - t134) * t173;
t53 = (t100 * t226 + t101 * t230) * t228 + t246;
t44 = -t155 * t289 + t171 * t290;
t43 = t153 * t289 - t171 * t291;
t40 = (t259 - t291) * t232 + t230 * t243;
t39 = t226 * t243 + t232 * t290 + t260;
t38 = t124 + t290 * t190 + (-t146 - t289) * t174;
t37 = t135 + t289 * t173 + (-t133 - t291) * t190;
t36 = -t153 * t290 + t155 * t291;
t33 = t232 * t87 + (t226 * t70 - t230 * t69) * t228;
t32 = t173 * t69 + t174 * t70 + t190 * t87;
t31 = (t226 * t291 + t230 * t290) * t228 + t246;
t26 = t123 + t291 * t174 + (-t134 - t290) * t173;
t25 = t232 * t75 + (t226 * t68 - t230 * t67) * t228;
t24 = t232 * t74 + (t226 * t66 - t230 * t65) * t228;
t23 = t173 * t67 + t174 * t68 + t190 * t75;
t22 = t173 * t65 + t174 * t66 + t190 * t74;
t102 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); t272 * t305 + t253 * t304 + (m(3) * t205 + m(4) * t167 - t307) * t283 + (m(3) * t204 + m(4) * t166 + t239) * t286 + t306 * (0.2e1 * t129 + 0.2e1 * t130 + t253); (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(7) + (t53 ^ 2 + t57 ^ 2 + t58 ^ 2) * m(6) + (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) * m(5) + m(4) * (t108 ^ 2 + t125 ^ 2 + t126 ^ 2) + m(3) * (t168 ^ 2 + t181 ^ 2 + t182 ^ 2) + (t8 + t18 + t25 + ((t161 * t212 + t163 * t193 + t165 * t194) * t226 - (t160 * t212 + t162 * t193 + t164 * t194) * t230) * t228 + (t199 * t286 + t201 * t223 + t203 * t224) * t286) * t286 + (-t7 - t17 - t24 - ((t161 * t211 + t163 * t191 + t165 * t192) * t226 - (t160 * t211 + t162 * t191 + t164 * t192) * t230) * t228 + (-t198 * t283 + t200 * t221 + t202 * t222) * t283 + (-t198 * t286 + t199 * t283 - t200 * t223 - t201 * t221 - t202 * t224 - t222 * t203) * t286) * t283 + (t12 + t21 + t33 + (t176 * t219 + t177 * t209 + t178 * t210 + t232 * t215) * t232 + (t176 * t212 + t177 * t193 + t178 * t194 + t215 * t286 + t216 * t223 + t217 * t224) * t286 + (-t176 * t211 - t177 * t191 - t178 * t192 + t215 * t283 - t216 * t221 - t217 * t222) * t283 + ((t161 * t219 + t163 * t209 + t165 * t210) * t226 - (t160 * t219 + t162 * t209 + t164 * t210) * t230 + (-t198 * t230 + t199 * t226 + t216 * t238 + t217 * t236) * t232 + ((t201 * t238 + t203 * t236) * t226 - (t200 * t238 + t202 * t236) * t230) * t228) * t228) * t232; t219 * t242; (m(4) * t108 + m(5) * t73 + m(6) * t53 + m(7) * t31) * t219 + (m(4) * t126 + m(5) * t85 + m(6) * t58 + m(7) * t40) * t212 + (m(4) * t125 + m(5) * t84 + m(6) * t57 + m(7) * t39) * t211; (t211 ^ 2 + t212 ^ 2 + t219 ^ 2) * t242; t239 * t174 + t307 * t173 + 0.2e1 * t306 * (-t173 * t134 + t123); (t26 * t31 + t37 * t40 + t38 * t39) * m(7) + (t53 * t54 + t57 * t62 + t58 * t61) * m(6) + (t73 * t86 + t84 * t89 + t85 * t88) * m(5) + (t32 / 0.2e1 + t262) * t232 + (t33 / 0.2e1 + t261) * t190 + (t25 / 0.2e1 + t264) * t174 + (t24 / 0.2e1 + t265) * t173 + ((-t22 / 0.2e1 - t267) * t230 + (t23 / 0.2e1 + t266) * t226) * t228; (m(5) * t86 + m(6) * t54 + m(7) * t26) * t219 + (m(5) * t88 + m(6) * t61 + m(7) * t37) * t212 + (m(5) * t89 + m(6) * t62 + m(7) * t38) * t211; (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) * m(6) + (t86 ^ 2 + t88 ^ 2 + t89 ^ 2) * m(5) + (t11 + t20 + t32) * t190 + (t6 + t16 + t23) * t174 + (t5 + t15 + t22) * t173; t153 * t240 + t155 * t241; (t31 * t36 + t39 * t44 + t40 * t43) * m(7) + (t53 * t64 + t57 * t72 + t58 * t71) * m(6) + t263 * t232 + t261 * t171 + t264 * t155 + t265 * t153 + (t226 * t268 - t230 * t269) * t228; (m(6) * t64 + m(7) * t36) * t219 + (m(6) * t71 + m(7) * t43) * t212 + (m(6) * t72 + m(7) * t44) * t211; (t26 * t36 + t37 * t43 + t38 * t44) * m(7) + (t54 * t64 + t61 * t71 + t62 * t72) * m(6) + t263 * t190 + t268 * t174 + t269 * t173 + t262 * t171 + t266 * t155 + t267 * t153; (t36 ^ 2 + t43 ^ 2 + t44 ^ 2) * m(7) + (t64 ^ 2 + t71 ^ 2 + t72 ^ 2) * m(6) + (t10 + t19) * t171 + (t4 + t14) * t155 + (t3 + t13) * t153; t46 * m(7); t7 * t298 + t8 * t297 + t232 * t299 + t12 * t296 + (t31 * t46 + t39 * t56 + t40 * t55) * m(7) + (-t230 * t1 / 0.2e1 + t226 * t300) * t228; (t211 * t56 + t212 * t55 + t219 * t46) * m(7); (t26 * t46 + t37 * t55 + t38 * t56) * m(7) + t11 * t296 + t190 * t299 + t6 * t297 + t174 * t300 + t173 * t301 + t5 * t298; t155 * t300 + t4 * t297 + t10 * t296 + t3 * t298 + t153 * t301 + t171 * t299 + (t36 * t46 + t43 * t55 + t44 * t56) * m(7); t138 * t1 + t149 * t9 + t140 * t2 + (t46 ^ 2 + t55 ^ 2 + t56 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t102(1) t102(2) t102(4) t102(7) t102(11) t102(16); t102(2) t102(3) t102(5) t102(8) t102(12) t102(17); t102(4) t102(5) t102(6) t102(9) t102(13) t102(18); t102(7) t102(8) t102(9) t102(10) t102(14) t102(19); t102(11) t102(12) t102(13) t102(14) t102(15) t102(20); t102(16) t102(17) t102(18) t102(19) t102(20) t102(21);];
Mq  = res;
