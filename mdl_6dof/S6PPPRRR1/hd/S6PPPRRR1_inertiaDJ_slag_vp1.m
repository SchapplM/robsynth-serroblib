% Calculate time derivative of joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPPRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:14
% EndTime: 2019-03-08 18:39:41
% DurationCPUTime: 17.94s
% Computational Cost: add. (189198->672), mult. (550712->964), div. (0->0), fcn. (719606->18), ass. (0->291)
t237 = sin(pkin(14));
t243 = cos(pkin(12));
t244 = cos(pkin(7));
t238 = sin(pkin(13));
t239 = sin(pkin(12));
t242 = cos(pkin(13));
t245 = cos(pkin(6));
t289 = t243 * t245;
t267 = -t238 * t239 + t242 * t289;
t257 = t267 * t244;
t265 = t238 * t289 + t239 * t242;
t240 = sin(pkin(7));
t241 = sin(pkin(6));
t293 = t240 * t241;
t296 = cos(pkin(14));
t226 = t265 * t296 + (-t243 * t293 + t257) * t237;
t291 = t241 * t244;
t233 = -t240 * t267 - t243 * t291;
t248 = sin(qJ(4));
t281 = t241 * t296;
t269 = t240 * t281;
t253 = t237 * t265 + t243 * t269 - t257 * t296;
t297 = cos(pkin(8));
t251 = t253 * t297;
t295 = sin(pkin(8));
t302 = cos(qJ(4));
t212 = t226 * t302 + (t233 * t295 - t251) * t248;
t221 = t233 * t297 + t253 * t295;
t247 = sin(qJ(5));
t301 = cos(qJ(5));
t193 = t212 * t301 + t221 * t247;
t268 = t302 * t295;
t211 = t226 * t248 - t233 * t268 + t251 * t302;
t205 = t211 * qJD(4);
t149 = qJD(5) * t193 - t205 * t247;
t263 = -t212 * t247 + t221 * t301;
t150 = qJD(5) * t263 - t205 * t301;
t206 = t212 * qJD(4);
t121 = rSges(6,1) * t150 - rSges(6,2) * t149 + rSges(6,3) * t206;
t294 = t239 * t245;
t266 = -t238 * t243 - t242 * t294;
t256 = t266 * t244;
t264 = -t238 * t294 + t242 * t243;
t227 = t264 * t296 + (t239 * t293 + t256) * t237;
t234 = t239 * t291 - t240 * t266;
t252 = t237 * t264 - t239 * t269 - t256 * t296;
t250 = t252 * t297;
t214 = t227 * t302 + (t234 * t295 - t250) * t248;
t222 = t234 * t297 + t252 * t295;
t195 = t214 * t301 + t222 * t247;
t213 = t227 * t248 - t234 * t268 + t250 * t302;
t207 = t213 * qJD(4);
t151 = qJD(5) * t195 - t207 * t247;
t262 = -t214 * t247 + t222 * t301;
t152 = qJD(5) * t262 - t207 * t301;
t208 = t214 * qJD(4);
t122 = rSges(6,1) * t152 - rSges(6,2) * t151 + rSges(6,3) * t208;
t138 = rSges(6,1) * t193 + rSges(6,2) * t263 + rSges(6,3) * t211;
t139 = rSges(6,1) * t195 + rSges(6,2) * t262 + rSges(6,3) * t213;
t57 = t121 * t213 - t122 * t211 + t138 * t208 - t139 * t206;
t320 = m(6) * t57;
t290 = t242 * t244;
t292 = t240 * t245;
t232 = t238 * t281 + (t241 * t290 + t292) * t237;
t236 = -t242 * t293 + t244 * t245;
t255 = t296 * t292 + (-t237 * t238 + t290 * t296) * t241;
t254 = t255 * t297;
t220 = t232 * t302 + (t236 * t295 + t254) * t248;
t216 = t220 * qJD(4);
t219 = t232 * t248 - t236 * t268 - t254 * t302;
t246 = sin(qJ(6));
t249 = cos(qJ(6));
t153 = -t193 * t246 + t211 * t249;
t154 = t193 * t249 + t211 * t246;
t103 = Icges(7,5) * t154 + Icges(7,6) * t153 - Icges(7,3) * t263;
t105 = Icges(7,4) * t154 + Icges(7,2) * t153 - Icges(7,6) * t263;
t107 = Icges(7,1) * t154 + Icges(7,4) * t153 - Icges(7,5) * t263;
t111 = -qJD(6) * t154 - t150 * t246 + t206 * t249;
t112 = qJD(6) * t153 + t150 * t249 + t206 * t246;
t78 = Icges(7,5) * t112 + Icges(7,6) * t111 + Icges(7,3) * t149;
t80 = Icges(7,4) * t112 + Icges(7,2) * t111 + Icges(7,6) * t149;
t82 = Icges(7,1) * t112 + Icges(7,4) * t111 + Icges(7,5) * t149;
t16 = t103 * t149 + t105 * t111 + t107 * t112 + t153 * t80 + t154 * t82 - t263 * t78;
t155 = -t195 * t246 + t213 * t249;
t156 = t195 * t249 + t213 * t246;
t104 = Icges(7,5) * t156 + Icges(7,6) * t155 - Icges(7,3) * t262;
t106 = Icges(7,4) * t156 + Icges(7,2) * t155 - Icges(7,6) * t262;
t108 = Icges(7,1) * t156 + Icges(7,4) * t155 - Icges(7,5) * t262;
t113 = -qJD(6) * t156 - t152 * t246 + t208 * t249;
t114 = qJD(6) * t155 + t152 * t249 + t208 * t246;
t79 = Icges(7,5) * t114 + Icges(7,6) * t113 + Icges(7,3) * t151;
t81 = Icges(7,4) * t114 + Icges(7,2) * t113 + Icges(7,6) * t151;
t83 = Icges(7,1) * t114 + Icges(7,4) * t113 + Icges(7,5) * t151;
t17 = t104 * t149 + t106 * t111 + t108 * t112 + t153 * t81 + t154 * t83 - t263 * t79;
t228 = t236 * t297 - t255 * t295;
t210 = t220 * t301 + t228 * t247;
t190 = -t210 * t246 + t219 * t249;
t191 = t210 * t249 + t219 * t246;
t261 = -t220 * t247 + t228 * t301;
t126 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t261;
t127 = Icges(7,4) * t191 + Icges(7,2) * t190 - Icges(7,6) * t261;
t128 = Icges(7,1) * t191 + Icges(7,4) * t190 - Icges(7,5) * t261;
t215 = t219 * qJD(4);
t189 = qJD(5) * t261 - t215 * t301;
t140 = -qJD(6) * t191 - t189 * t246 + t216 * t249;
t141 = qJD(6) * t190 + t189 * t249 + t216 * t246;
t188 = qJD(5) * t210 - t215 * t247;
t97 = Icges(7,5) * t141 + Icges(7,6) * t140 + Icges(7,3) * t188;
t98 = Icges(7,4) * t141 + Icges(7,2) * t140 + Icges(7,6) * t188;
t99 = Icges(7,1) * t141 + Icges(7,4) * t140 + Icges(7,5) * t188;
t32 = t111 * t127 + t112 * t128 + t126 * t149 + t153 * t98 + t154 * t99 - t263 * t97;
t52 = -t103 * t263 + t105 * t153 + t107 * t154;
t53 = -t104 * t263 + t106 * t153 + t108 * t154;
t65 = -t126 * t263 + t127 * t153 + t128 * t154;
t3 = t16 * t211 + t17 * t213 + t206 * t52 + t208 * t53 + t216 * t65 + t219 * t32;
t115 = Icges(6,5) * t150 - Icges(6,6) * t149 + Icges(6,3) * t206;
t117 = Icges(6,4) * t150 - Icges(6,2) * t149 + Icges(6,6) * t206;
t119 = Icges(6,1) * t150 - Icges(6,4) * t149 + Icges(6,5) * t206;
t132 = Icges(6,5) * t193 + Icges(6,6) * t263 + Icges(6,3) * t211;
t134 = Icges(6,4) * t193 + Icges(6,2) * t263 + Icges(6,6) * t211;
t136 = Icges(6,1) * t193 + Icges(6,4) * t263 + Icges(6,5) * t211;
t41 = t115 * t211 + t117 * t263 + t119 * t193 + t132 * t206 - t134 * t149 + t136 * t150;
t116 = Icges(6,5) * t152 - Icges(6,6) * t151 + Icges(6,3) * t208;
t118 = Icges(6,4) * t152 - Icges(6,2) * t151 + Icges(6,6) * t208;
t120 = Icges(6,1) * t152 - Icges(6,4) * t151 + Icges(6,5) * t208;
t133 = Icges(6,5) * t195 + Icges(6,6) * t262 + Icges(6,3) * t213;
t135 = Icges(6,4) * t195 + Icges(6,2) * t262 + Icges(6,6) * t213;
t137 = Icges(6,1) * t195 + Icges(6,4) * t262 + Icges(6,5) * t213;
t42 = t116 * t211 + t118 * t263 + t120 * t193 + t133 * t206 - t135 * t149 + t137 * t150;
t142 = Icges(6,5) * t189 - Icges(6,6) * t188 + Icges(6,3) * t216;
t143 = Icges(6,4) * t189 - Icges(6,2) * t188 + Icges(6,6) * t216;
t144 = Icges(6,1) * t189 - Icges(6,4) * t188 + Icges(6,5) * t216;
t160 = Icges(6,5) * t210 + Icges(6,6) * t261 + Icges(6,3) * t219;
t161 = Icges(6,4) * t210 + Icges(6,2) * t261 + Icges(6,6) * t219;
t162 = Icges(6,1) * t210 + Icges(6,4) * t261 + Icges(6,5) * t219;
t49 = t142 * t211 + t143 * t263 + t144 * t193 - t149 * t161 + t150 * t162 + t160 * t206;
t72 = t132 * t211 + t134 * t263 + t136 * t193;
t73 = t133 * t211 + t135 * t263 + t137 * t193;
t91 = t160 * t211 + t161 * t263 + t162 * t193;
t319 = t206 * t72 + t208 * t73 + t211 * t41 + t213 * t42 + t216 * t91 + t219 * t49 + t3;
t18 = t103 * t151 + t105 * t113 + t107 * t114 + t155 * t80 + t156 * t82 - t262 * t78;
t19 = t104 * t151 + t106 * t113 + t108 * t114 + t155 * t81 + t156 * t83 - t262 * t79;
t33 = t113 * t127 + t114 * t128 + t126 * t151 + t155 * t98 + t156 * t99 - t262 * t97;
t54 = -t103 * t262 + t105 * t155 + t107 * t156;
t55 = -t104 * t262 + t106 * t155 + t108 * t156;
t66 = -t126 * t262 + t127 * t155 + t128 * t156;
t4 = t18 * t211 + t19 * t213 + t206 * t54 + t208 * t55 + t216 * t66 + t219 * t33;
t43 = t115 * t213 + t117 * t262 + t119 * t195 + t132 * t208 - t134 * t151 + t136 * t152;
t44 = t116 * t213 + t118 * t262 + t120 * t195 + t133 * t208 - t135 * t151 + t137 * t152;
t50 = t142 * t213 + t143 * t262 + t144 * t195 - t151 * t161 + t152 * t162 + t160 * t208;
t74 = t132 * t213 + t134 * t262 + t136 * t195;
t75 = t133 * t213 + t135 * t262 + t137 * t195;
t92 = t160 * t213 + t161 * t262 + t162 * t195;
t318 = t206 * t74 + t208 * t75 + t211 * t43 + t213 * t44 + t216 * t92 + t219 * t50 + t4;
t96 = t160 * t219 + t161 * t261 + t162 * t210;
t300 = t216 * t96;
t47 = t115 * t219 + t117 * t261 + t119 * t210 + t132 * t216 - t134 * t188 + t136 * t189;
t48 = t116 * t219 + t118 * t261 + t120 * t210 + t133 * t216 - t135 * t188 + t137 * t189;
t56 = t142 * t219 + t143 * t261 + t144 * t210 + t160 * t216 - t161 * t188 + t162 * t189;
t20 = t103 * t188 + t105 * t140 + t107 * t141 + t190 * t80 + t191 * t82 - t261 * t78;
t21 = t104 * t188 + t106 * t140 + t108 * t141 + t190 * t81 + t191 * t83 - t261 * t79;
t36 = t126 * t188 + t127 * t140 + t128 * t141 + t190 * t98 + t191 * t99 - t261 * t97;
t58 = -t103 * t261 + t105 * t190 + t107 * t191;
t59 = -t104 * t261 + t106 * t190 + t108 * t191;
t70 = -t126 * t261 + t127 * t190 + t128 * t191;
t6 = t20 * t211 + t206 * t58 + t208 * t59 + t21 * t213 + t216 * t70 + t219 * t36;
t86 = t132 * t219 + t134 * t261 + t136 * t210;
t87 = t133 * t219 + t135 * t261 + t137 * t210;
t317 = t206 * t86 + t208 * t87 + t211 * t47 + t213 * t48 + t219 * t56 + t300 + t6;
t110 = rSges(7,1) * t156 + rSges(7,2) * t155 - rSges(7,3) * t262;
t286 = pkin(5) * t195 - pkin(11) * t262 + t110;
t109 = rSges(7,1) * t154 + rSges(7,2) * t153 - rSges(7,3) * t263;
t287 = pkin(5) * t193 - pkin(11) * t263 + t109;
t316 = -t206 * t286 + t208 * t287;
t315 = m(7) / 0.2e1;
t314 = t149 / 0.2e1;
t313 = t151 / 0.2e1;
t312 = t188 / 0.2e1;
t311 = -t263 / 0.2e1;
t310 = -t262 / 0.2e1;
t309 = t206 / 0.2e1;
t308 = t208 / 0.2e1;
t307 = -t261 / 0.2e1;
t306 = t216 / 0.2e1;
t305 = t221 / 0.2e1;
t304 = t222 / 0.2e1;
t303 = t228 / 0.2e1;
t84 = rSges(7,1) * t112 + rSges(7,2) * t111 + rSges(7,3) * t149;
t299 = pkin(5) * t150 + pkin(11) * t149 + t84;
t85 = rSges(7,1) * t114 + rSges(7,2) * t113 + rSges(7,3) * t151;
t298 = pkin(5) * t152 + pkin(11) * t151 + t85;
t100 = rSges(7,1) * t141 + rSges(7,2) * t140 + rSges(7,3) * t188;
t288 = pkin(5) * t189 + pkin(11) * t188 + t100;
t129 = rSges(7,1) * t191 + rSges(7,2) * t190 - rSges(7,3) * t261;
t285 = pkin(5) * t210 - pkin(11) * t261 + t129;
t279 = 0.2e1 * m(6);
t277 = 0.2e1 * m(7);
t274 = 0.2e1 * t299;
t273 = 0.2e1 * t298;
t29 = -t211 * t298 + t213 * t299 + t316;
t272 = m(7) * t29 + t320;
t34 = t206 * t285 + t211 * t288 - t216 * t287 - t219 * t299;
t145 = rSges(6,1) * t189 - rSges(6,2) * t188 + rSges(6,3) * t216;
t163 = rSges(6,1) * t210 + rSges(6,2) * t261 + rSges(6,3) * t219;
t63 = -t121 * t219 - t138 * t216 + t145 * t211 + t163 * t206;
t271 = m(6) * t63 + m(7) * t34;
t35 = -t208 * t285 - t213 * t288 + t216 * t286 + t219 * t298;
t64 = t122 * t219 + t139 * t216 - t145 * t213 - t163 * t208;
t270 = m(6) * t64 + m(7) * t35;
t179 = -rSges(5,1) * t205 - rSges(5,2) * t206;
t180 = -rSges(5,1) * t207 - rSges(5,2) * t208;
t125 = t179 * t222 - t180 * t221;
t181 = -pkin(4) * t205 + pkin(10) * t206;
t159 = t222 * t181;
t182 = -pkin(4) * t207 + pkin(10) * t208;
t38 = t159 + t299 * t222 + (-t182 - t298) * t221;
t69 = t121 * t222 + t159 + (-t122 - t182) * t221;
t260 = m(5) * t125 + m(6) * t69 + m(7) * t38;
t202 = -rSges(5,1) * t215 - rSges(5,2) * t216;
t130 = -t179 * t228 + t202 * t221;
t203 = -pkin(4) * t215 + pkin(10) * t216;
t186 = t221 * t203;
t45 = t186 + t288 * t221 + (-t181 - t299) * t228;
t76 = t145 * t221 + t186 + (-t121 - t181) * t228;
t259 = m(5) * t130 + m(6) * t76 + m(7) * t45;
t131 = t180 * t228 - t202 * t222;
t171 = t228 * t182;
t46 = t171 + t298 * t228 + (-t203 - t288) * t222;
t77 = t122 * t228 + t171 + (-t145 - t203) * t222;
t258 = m(5) * t131 + m(6) * t77 + m(7) * t46;
t37 = t109 * t151 - t110 * t149 - t262 * t84 + t263 * t85;
t204 = pkin(4) * t220 + pkin(10) * t219;
t201 = -Icges(5,1) * t215 - Icges(5,4) * t216;
t200 = -Icges(5,4) * t215 - Icges(5,2) * t216;
t199 = -Icges(5,5) * t215 - Icges(5,6) * t216;
t198 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t228;
t197 = Icges(5,1) * t220 - Icges(5,4) * t219 + Icges(5,5) * t228;
t196 = Icges(5,4) * t220 - Icges(5,2) * t219 + Icges(5,6) * t228;
t187 = t221 * t204;
t185 = pkin(4) * t214 + pkin(10) * t213;
t184 = pkin(4) * t212 + pkin(10) * t211;
t178 = -Icges(5,1) * t207 - Icges(5,4) * t208;
t177 = -Icges(5,1) * t205 - Icges(5,4) * t206;
t176 = -Icges(5,4) * t207 - Icges(5,2) * t208;
t175 = -Icges(5,4) * t205 - Icges(5,2) * t206;
t174 = -Icges(5,5) * t207 - Icges(5,6) * t208;
t173 = -Icges(5,5) * t205 - Icges(5,6) * t206;
t172 = t228 * t185;
t170 = t222 * t184;
t169 = rSges(5,1) * t214 - rSges(5,2) * t213 + rSges(5,3) * t222;
t168 = rSges(5,1) * t212 - rSges(5,2) * t211 + rSges(5,3) * t221;
t167 = Icges(5,1) * t214 - Icges(5,4) * t213 + Icges(5,5) * t222;
t166 = Icges(5,1) * t212 - Icges(5,4) * t211 + Icges(5,5) * t221;
t165 = Icges(5,4) * t214 - Icges(5,2) * t213 + Icges(5,6) * t222;
t164 = Icges(5,4) * t212 - Icges(5,2) * t211 + Icges(5,6) * t221;
t102 = t139 * t219 - t163 * t213;
t101 = -t138 * t219 + t163 * t211;
t95 = t138 * t213 - t139 * t211;
t94 = t139 * t228 + t172 + (-t163 - t204) * t222;
t93 = t163 * t221 + t187 + (-t138 - t184) * t228;
t90 = -t110 * t261 + t129 * t262;
t89 = t109 * t261 - t129 * t263;
t88 = t138 * t222 + t170 + (-t139 - t185) * t221;
t71 = -t109 * t262 + t110 * t263;
t68 = -t213 * t285 + t219 * t286;
t67 = t211 * t285 - t219 * t287;
t62 = t172 + t286 * t228 + (-t204 - t285) * t222;
t61 = t187 + t285 * t221 + (-t184 - t287) * t228;
t60 = -t211 * t286 + t213 * t287;
t51 = t170 + t287 * t222 + (-t185 - t286) * t221;
t40 = t100 * t262 + t110 * t188 - t129 * t151 - t261 * t85;
t39 = -t100 * t263 - t109 * t188 + t129 * t149 + t261 * t84;
t31 = t221 * t58 + t222 * t59 + t228 * t70;
t30 = t211 * t58 + t213 * t59 + t219 * t70;
t28 = -t261 * t70 - t262 * t59 - t263 * t58;
t27 = t221 * t54 + t222 * t55 + t228 * t66;
t26 = t221 * t52 + t222 * t53 + t228 * t65;
t25 = t211 * t54 + t213 * t55 + t219 * t66;
t24 = t211 * t52 + t213 * t53 + t219 * t65;
t23 = -t261 * t66 - t262 * t55 - t263 * t54;
t22 = -t261 * t65 - t262 * t53 - t263 * t52;
t15 = t221 * t47 + t222 * t48 + t228 * t56;
t14 = t221 * t43 + t222 * t44 + t228 * t50;
t13 = t221 * t41 + t222 * t42 + t228 * t49;
t9 = t20 * t221 + t21 * t222 + t228 * t36;
t8 = t18 * t221 + t19 * t222 + t228 * t33;
t7 = t16 * t221 + t17 * t222 + t228 * t32;
t5 = t149 * t58 + t151 * t59 + t188 * t70 - t20 * t263 - t21 * t262 - t261 * t36;
t2 = t149 * t54 + t151 * t55 - t18 * t263 + t188 * t66 - t19 * t262 - t261 * t33;
t1 = t149 * t52 + t151 * t53 - t16 * t263 - t17 * t262 + t188 * t65 - t261 * t32;
t10 = [0; 0; 0; 0; 0; 0; (m(5) * t179 + m(6) * t121 + t274 * t315) * t222 + (-m(5) * t180 - m(6) * t122 - t273 * t315) * t221 + 0.2e1 * (m(6) / 0.2e1 + t315) * (-t221 * t182 + t159); t260 * t245 + (t239 * t259 - t243 * t258) * t241; t233 * t258 + t234 * t259 + t236 * t260; (t38 * t51 + t45 * t61 + t46 * t62) * t277 + t221 * t7 + t228 * t9 + t222 * t8 + t228 * t15 + t221 * t13 + t222 * t14 + (t69 * t88 + t76 * t93 + t77 * t94) * t279 + t222 * ((-t165 * t208 - t167 * t207 + t174 * t222 - t176 * t213 + t178 * t214) * t222 + (-t164 * t208 - t166 * t207 + t173 * t222 - t175 * t213 + t177 * t214) * t221 + (-t196 * t208 - t197 * t207 + t199 * t222 - t200 * t213 + t201 * t214) * t228) + t228 * ((-t165 * t216 - t167 * t215 + t174 * t228 - t176 * t219 + t178 * t220) * t222 + (-t164 * t216 - t166 * t215 + t173 * t228 - t175 * t219 + t177 * t220) * t221 + (-t196 * t216 - t197 * t215 + t199 * t228 - t200 * t219 + t201 * t220) * t228) + t221 * ((-t165 * t206 - t167 * t205 + t174 * t221 - t176 * t211 + t178 * t212) * t222 + (-t164 * t206 - t166 * t205 + t173 * t221 - t175 * t211 + t177 * t212) * t221 + (-t196 * t206 - t197 * t205 + t199 * t221 - t200 * t211 + t201 * t212) * t228) + 0.2e1 * m(5) * ((-t168 * t228 + t198 * t221) * t130 + (t169 * t228 - t198 * t222) * t131 + (t168 * t222 - t169 * t221) * t125); t320 + (-t211 * t273 + t213 * t274 + 0.2e1 * t316) * t315; t272 * t245 + (t239 * t271 - t243 * t270) * t241; t233 * t270 + t234 * t271 + t236 * t272; (t29 * t51 + t34 * t61 + t35 * t62 + t38 * t60 + t45 * t67 + t46 * t68) * m(7) + (t101 * t76 + t102 * t77 + t57 * t88 + t63 * t93 + t64 * t94 + t69 * t95) * m(6) + (t9 / 0.2e1 + t15 / 0.2e1) * t219 + (t8 / 0.2e1 + t14 / 0.2e1) * t213 + (t7 / 0.2e1 + t13 / 0.2e1) * t211 + (t221 * t72 + t222 * t73 + t228 * t91 + t26) * t309 + (t221 * t74 + t222 * t75 + t228 * t92 + t27) * t308 + (t221 * t86 + t222 * t87 + t228 * t96 + t31) * t306 + t319 * t305 + t318 * t304 + t317 * t303; (t29 * t60 + t34 * t67 + t35 * t68) * t277 + t216 * t30 + (t101 * t63 + t102 * t64 + t57 * t95) * t279 + (t300 + t317) * t219 + (t216 * t87 + t318) * t213 + (t216 * t86 + t319) * t211 + (t211 * t74 + t213 * t75 + t219 * t92 + t25) * t208 + (t211 * t72 + t213 * t73 + t219 * t91 + t24) * t206; t37 * m(7); 0.2e1 * (t245 * t37 + (t239 * t39 - t243 * t40) * t241) * t315; (t233 * t40 + t234 * t39 + t236 * t37) * m(7); (t37 * t51 + t38 * t71 + t39 * t61 + t40 * t62 + t45 * t89 + t46 * t90) * m(7) + t1 * t305 + t5 * t303 + t2 * t304 + t31 * t312 + t9 * t307 + t26 * t314 + t7 * t311 + t27 * t313 + t8 * t310; t25 * t313 + t4 * t310 + t23 * t308 + t213 * t2 / 0.2e1 + t24 * t314 + t3 * t311 + (t29 * t71 + t34 * t89 + t35 * t90 + t37 * t60 + t39 * t67 + t40 * t68) * m(7) + t30 * t312 + t6 * t307 + t22 * t309 + t211 * t1 / 0.2e1 + t28 * t306 + t219 * t5 / 0.2e1; t149 * t22 - t263 * t1 + t188 * t28 - t261 * t5 + t151 * t23 - t262 * t2 + (t37 * t71 + t39 * t89 + t40 * t90) * t277;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
