% Calculate time derivative of joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:40
% EndTime: 2019-12-05 15:19:01
% DurationCPUTime: 10.00s
% Computational Cost: add. (68194->651), mult. (200535->926), div. (0->0), fcn. (249586->14), ass. (0->268)
t234 = sin(pkin(10));
t237 = cos(pkin(5));
t233 = sin(pkin(11));
t236 = cos(pkin(10));
t274 = t236 * t233;
t278 = cos(pkin(11));
t227 = t234 * t278 + t237 * t274;
t240 = sin(qJ(3));
t263 = t237 * t278;
t275 = t234 * t233;
t245 = -t236 * t263 + t275;
t279 = cos(pkin(6));
t243 = t245 * t279;
t235 = sin(pkin(5));
t277 = sin(pkin(6));
t264 = t235 * t277;
t284 = cos(qJ(3));
t210 = t227 * t284 + (-t236 * t264 - t243) * t240;
t265 = t235 * t279;
t221 = -t236 * t265 + t245 * t277;
t239 = sin(qJ(4));
t283 = cos(qJ(4));
t191 = t210 * t283 + t221 * t239;
t253 = t284 * t277;
t251 = t235 * t253;
t209 = t227 * t240 + t236 * t251 + t284 * t243;
t205 = t209 * qJD(3);
t149 = qJD(4) * t191 - t205 * t239;
t249 = -t210 * t239 + t221 * t283;
t150 = t249 * qJD(4) - t205 * t283;
t206 = t210 * qJD(3);
t121 = rSges(5,1) * t150 - rSges(5,2) * t149 + rSges(5,3) * t206;
t228 = t236 * t278 - t237 * t275;
t244 = t234 * t263 + t274;
t242 = t244 * t279;
t212 = t228 * t284 + (t234 * t264 - t242) * t240;
t222 = t234 * t265 + t244 * t277;
t193 = t212 * t283 + t222 * t239;
t211 = t228 * t240 - t234 * t251 + t284 * t242;
t207 = t211 * qJD(3);
t151 = qJD(4) * t193 - t207 * t239;
t248 = -t212 * t239 + t222 * t283;
t152 = t248 * qJD(4) - t207 * t283;
t208 = t212 * qJD(3);
t122 = rSges(5,1) * t152 - rSges(5,2) * t151 + rSges(5,3) * t208;
t132 = rSges(5,1) * t191 + rSges(5,2) * t249 + rSges(5,3) * t209;
t133 = rSges(5,1) * t193 + rSges(5,2) * t248 + rSges(5,3) * t211;
t57 = t121 * t211 - t122 * t209 + t132 * t208 - t133 * t206;
t304 = m(5) * t57;
t252 = t279 * t278;
t303 = (-t233 * t240 + t284 * t252) * t235 + t237 * t253;
t220 = t237 * t277 * t240 + (t284 * t233 + t240 * t252) * t235;
t216 = t220 * qJD(3);
t238 = sin(qJ(5));
t241 = cos(qJ(5));
t153 = -t191 * t238 + t209 * t241;
t154 = t191 * t241 + t209 * t238;
t103 = Icges(6,5) * t154 + Icges(6,6) * t153 - Icges(6,3) * t249;
t105 = Icges(6,4) * t154 + Icges(6,2) * t153 - Icges(6,6) * t249;
t107 = Icges(6,1) * t154 + Icges(6,4) * t153 - Icges(6,5) * t249;
t111 = -qJD(5) * t154 - t150 * t238 + t206 * t241;
t112 = qJD(5) * t153 + t150 * t241 + t206 * t238;
t76 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t149;
t78 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t149;
t80 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t149;
t16 = t103 * t149 + t105 * t111 + t107 * t112 + t153 * t78 + t154 * t80 - t249 * t76;
t155 = -t193 * t238 + t211 * t241;
t156 = t193 * t241 + t211 * t238;
t104 = Icges(6,5) * t156 + Icges(6,6) * t155 - Icges(6,3) * t248;
t106 = Icges(6,4) * t156 + Icges(6,2) * t155 - Icges(6,6) * t248;
t108 = Icges(6,1) * t156 + Icges(6,4) * t155 - Icges(6,5) * t248;
t113 = -qJD(5) * t156 - t152 * t238 + t208 * t241;
t114 = qJD(5) * t155 + t152 * t241 + t208 * t238;
t77 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t151;
t79 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t151;
t81 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t151;
t17 = t104 * t149 + t106 * t111 + t108 * t112 + t153 * t79 + t154 * t81 - t249 * t77;
t226 = t237 * t279 - t278 * t264;
t214 = t220 * t283 + t226 * t239;
t194 = -t214 * t238 - t241 * t303;
t195 = t214 * t241 - t238 * t303;
t247 = -t220 * t239 + t226 * t283;
t134 = Icges(6,5) * t195 + Icges(6,6) * t194 - Icges(6,3) * t247;
t135 = Icges(6,4) * t195 + Icges(6,2) * t194 - Icges(6,6) * t247;
t136 = Icges(6,1) * t195 + Icges(6,4) * t194 - Icges(6,5) * t247;
t215 = t303 * qJD(3);
t189 = t247 * qJD(4) + t215 * t283;
t140 = -qJD(5) * t195 - t189 * t238 + t216 * t241;
t141 = qJD(5) * t194 + t189 * t241 + t216 * t238;
t188 = qJD(4) * t214 + t215 * t239;
t97 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t188;
t98 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t188;
t99 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t188;
t32 = t111 * t135 + t112 * t136 + t134 * t149 + t153 * t98 + t154 * t99 - t249 * t97;
t51 = -t103 * t249 + t105 * t153 + t107 * t154;
t52 = -t104 * t249 + t106 * t153 + t108 * t154;
t65 = -t134 * t249 + t135 * t153 + t136 * t154;
t3 = t16 * t209 + t17 * t211 + t206 * t51 + t208 * t52 + t216 * t65 - t303 * t32;
t115 = Icges(5,5) * t150 - Icges(5,6) * t149 + Icges(5,3) * t206;
t117 = Icges(5,4) * t150 - Icges(5,2) * t149 + Icges(5,6) * t206;
t119 = Icges(5,1) * t150 - Icges(5,4) * t149 + Icges(5,5) * t206;
t126 = Icges(5,5) * t191 + Icges(5,6) * t249 + Icges(5,3) * t209;
t128 = Icges(5,4) * t191 + Icges(5,2) * t249 + Icges(5,6) * t209;
t130 = Icges(5,1) * t191 + Icges(5,4) * t249 + Icges(5,5) * t209;
t41 = t115 * t209 + t117 * t249 + t119 * t191 + t126 * t206 - t128 * t149 + t130 * t150;
t116 = Icges(5,5) * t152 - Icges(5,6) * t151 + Icges(5,3) * t208;
t118 = Icges(5,4) * t152 - Icges(5,2) * t151 + Icges(5,6) * t208;
t120 = Icges(5,1) * t152 - Icges(5,4) * t151 + Icges(5,5) * t208;
t127 = Icges(5,5) * t193 + Icges(5,6) * t248 + Icges(5,3) * t211;
t129 = Icges(5,4) * t193 + Icges(5,2) * t248 + Icges(5,6) * t211;
t131 = Icges(5,1) * t193 + Icges(5,4) * t248 + Icges(5,5) * t211;
t42 = t116 * t209 + t118 * t249 + t120 * t191 + t127 * t206 - t129 * t149 + t131 * t150;
t142 = Icges(5,5) * t189 - Icges(5,6) * t188 + Icges(5,3) * t216;
t143 = Icges(5,4) * t189 - Icges(5,2) * t188 + Icges(5,6) * t216;
t144 = Icges(5,1) * t189 - Icges(5,4) * t188 + Icges(5,5) * t216;
t166 = Icges(5,5) * t214 + Icges(5,6) * t247 - Icges(5,3) * t303;
t167 = Icges(5,4) * t214 + Icges(5,2) * t247 - Icges(5,6) * t303;
t168 = Icges(5,1) * t214 + Icges(5,4) * t247 - Icges(5,5) * t303;
t49 = t142 * t209 + t143 * t249 + t144 * t191 - t149 * t167 + t150 * t168 + t166 * t206;
t71 = t126 * t209 + t128 * t249 + t130 * t191;
t72 = t127 * t209 + t129 * t249 + t131 * t191;
t91 = t166 * t209 + t167 * t249 + t168 * t191;
t302 = t206 * t71 + t208 * t72 + t209 * t41 + t211 * t42 + t216 * t91 - t303 * t49 + t3;
t18 = t103 * t151 + t105 * t113 + t107 * t114 + t155 * t78 + t156 * t80 - t248 * t76;
t19 = t104 * t151 + t106 * t113 + t108 * t114 + t155 * t79 + t156 * t81 - t248 * t77;
t33 = t113 * t135 + t114 * t136 + t134 * t151 + t155 * t98 + t156 * t99 - t248 * t97;
t53 = -t103 * t248 + t105 * t155 + t107 * t156;
t54 = -t104 * t248 + t106 * t155 + t108 * t156;
t66 = -t134 * t248 + t135 * t155 + t136 * t156;
t4 = t18 * t209 + t19 * t211 + t206 * t53 + t208 * t54 + t216 * t66 - t303 * t33;
t43 = t115 * t211 + t117 * t248 + t119 * t193 + t126 * t208 - t128 * t151 + t130 * t152;
t44 = t116 * t211 + t118 * t248 + t120 * t193 + t127 * t208 - t129 * t151 + t131 * t152;
t50 = t142 * t211 + t143 * t248 + t144 * t193 - t151 * t167 + t152 * t168 + t166 * t208;
t73 = t126 * t211 + t128 * t248 + t130 * t193;
t74 = t127 * t211 + t129 * t248 + t131 * t193;
t92 = t166 * t211 + t167 * t248 + t168 * t193;
t301 = t206 * t73 + t208 * t74 + t209 * t43 + t211 * t44 + t216 * t92 - t303 * t50 + t4;
t96 = -t166 * t303 + t167 * t247 + t168 * t214;
t282 = t216 * t96;
t47 = -t115 * t303 + t117 * t247 + t119 * t214 + t126 * t216 - t128 * t188 + t130 * t189;
t48 = -t116 * t303 + t118 * t247 + t120 * t214 + t127 * t216 - t129 * t188 + t131 * t189;
t56 = -t142 * t303 + t143 * t247 + t144 * t214 + t166 * t216 - t167 * t188 + t168 * t189;
t20 = t103 * t188 + t105 * t140 + t107 * t141 + t194 * t78 + t195 * t80 - t247 * t76;
t21 = t104 * t188 + t106 * t140 + t108 * t141 + t194 * t79 + t195 * t81 - t247 * t77;
t36 = t134 * t188 + t135 * t140 + t136 * t141 + t194 * t98 + t195 * t99 - t247 * t97;
t58 = -t103 * t247 + t105 * t194 + t107 * t195;
t59 = -t104 * t247 + t106 * t194 + t108 * t195;
t75 = -t134 * t247 + t135 * t194 + t136 * t195;
t6 = t20 * t209 + t206 * t58 + t208 * t59 + t21 * t211 + t216 * t75 - t303 * t36;
t86 = -t126 * t303 + t128 * t247 + t130 * t214;
t87 = -t127 * t303 + t129 * t247 + t131 * t214;
t300 = t206 * t86 + t208 * t87 + t209 * t47 + t211 * t48 - t303 * t56 + t282 + t6;
t110 = rSges(6,1) * t156 + rSges(6,2) * t155 - rSges(6,3) * t248;
t271 = pkin(4) * t193 - pkin(9) * t248 + t110;
t109 = rSges(6,1) * t154 + rSges(6,2) * t153 - rSges(6,3) * t249;
t272 = pkin(4) * t191 - pkin(9) * t249 + t109;
t298 = -t271 * t206 + t272 * t208;
t297 = m(6) / 0.2e1;
t296 = t149 / 0.2e1;
t295 = t151 / 0.2e1;
t294 = t188 / 0.2e1;
t293 = -t249 / 0.2e1;
t292 = -t248 / 0.2e1;
t291 = t206 / 0.2e1;
t290 = t208 / 0.2e1;
t289 = -t247 / 0.2e1;
t288 = t216 / 0.2e1;
t287 = t221 / 0.2e1;
t286 = t222 / 0.2e1;
t285 = t226 / 0.2e1;
t82 = rSges(6,1) * t112 + rSges(6,2) * t111 + rSges(6,3) * t149;
t281 = pkin(4) * t150 + pkin(9) * t149 + t82;
t83 = rSges(6,1) * t114 + rSges(6,2) * t113 + rSges(6,3) * t151;
t280 = pkin(4) * t152 + pkin(9) * t151 + t83;
t100 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t188;
t273 = pkin(4) * t189 + pkin(9) * t188 + t100;
t139 = rSges(6,1) * t195 + rSges(6,2) * t194 - rSges(6,3) * t247;
t270 = pkin(4) * t214 - pkin(9) * t247 + t139;
t261 = 0.2e1 * m(5);
t259 = 0.2e1 * m(6);
t256 = 0.2e1 * t281;
t255 = 0.2e1 * t280;
t37 = t109 * t151 - t110 * t149 - t248 * t82 + t249 * t83;
t204 = pkin(3) * t220 - pkin(8) * t303;
t203 = pkin(3) * t215 + pkin(8) * t216;
t202 = rSges(4,1) * t215 - rSges(4,2) * t216;
t201 = Icges(4,1) * t215 - Icges(4,4) * t216;
t200 = Icges(4,4) * t215 - Icges(4,2) * t216;
t199 = Icges(4,5) * t215 - Icges(4,6) * t216;
t198 = rSges(4,1) * t220 + rSges(4,2) * t303 + rSges(4,3) * t226;
t197 = Icges(4,1) * t220 + Icges(4,4) * t303 + Icges(4,5) * t226;
t196 = Icges(4,4) * t220 + Icges(4,2) * t303 + Icges(4,6) * t226;
t187 = t221 * t204;
t186 = t221 * t203;
t184 = pkin(3) * t212 + pkin(8) * t211;
t183 = pkin(3) * t210 + pkin(8) * t209;
t182 = -pkin(3) * t207 + pkin(8) * t208;
t181 = -pkin(3) * t205 + pkin(8) * t206;
t180 = -rSges(4,1) * t207 - rSges(4,2) * t208;
t179 = -rSges(4,1) * t205 - rSges(4,2) * t206;
t178 = -Icges(4,1) * t207 - Icges(4,4) * t208;
t177 = -Icges(4,1) * t205 - Icges(4,4) * t206;
t176 = -Icges(4,4) * t207 - Icges(4,2) * t208;
t175 = -Icges(4,4) * t205 - Icges(4,2) * t206;
t174 = -Icges(4,5) * t207 - Icges(4,6) * t208;
t173 = -Icges(4,5) * t205 - Icges(4,6) * t206;
t172 = t226 * t184;
t171 = t226 * t182;
t170 = t222 * t183;
t169 = rSges(5,1) * t214 + rSges(5,2) * t247 - rSges(5,3) * t303;
t165 = rSges(4,1) * t212 - rSges(4,2) * t211 + rSges(4,3) * t222;
t164 = rSges(4,1) * t210 - rSges(4,2) * t209 + rSges(4,3) * t221;
t163 = Icges(4,1) * t212 - Icges(4,4) * t211 + Icges(4,5) * t222;
t162 = Icges(4,1) * t210 - Icges(4,4) * t209 + Icges(4,5) * t221;
t161 = Icges(4,4) * t212 - Icges(4,2) * t211 + Icges(4,6) * t222;
t160 = Icges(4,4) * t210 - Icges(4,2) * t209 + Icges(4,6) * t221;
t159 = t222 * t181;
t145 = rSges(5,1) * t189 - rSges(5,2) * t188 + rSges(5,3) * t216;
t138 = t180 * t226 - t202 * t222;
t137 = -t179 * t226 + t202 * t221;
t125 = t179 * t222 - t180 * t221;
t102 = -t133 * t303 - t169 * t211;
t101 = t132 * t303 + t169 * t209;
t95 = t132 * t211 - t133 * t209;
t94 = t133 * t226 + t172 + (-t169 - t204) * t222;
t93 = t169 * t221 + t187 + (-t132 - t183) * t226;
t90 = -t110 * t247 + t139 * t248;
t89 = t109 * t247 - t139 * t249;
t88 = t132 * t222 + t170 + (-t133 - t184) * t221;
t85 = t122 * t226 + t171 + (-t145 - t203) * t222;
t84 = t145 * t221 + t186 + (-t121 - t181) * t226;
t70 = -t109 * t248 + t110 * t249;
t69 = -t270 * t211 - t271 * t303;
t68 = t270 * t209 + t272 * t303;
t67 = t121 * t222 + t159 + (-t122 - t182) * t221;
t64 = -t122 * t303 + t133 * t216 - t145 * t211 - t169 * t208;
t63 = t121 * t303 - t132 * t216 + t145 * t209 + t169 * t206;
t62 = t172 + t271 * t226 + (-t204 - t270) * t222;
t61 = t187 + t270 * t221 + (-t183 - t272) * t226;
t60 = -t271 * t209 + t272 * t211;
t55 = t170 + t272 * t222 + (-t184 - t271) * t221;
t46 = t171 + t280 * t226 + (-t203 - t273) * t222;
t45 = t186 + t273 * t221 + (-t181 - t281) * t226;
t40 = t100 * t248 + t110 * t188 - t139 * t151 - t247 * t83;
t39 = -t100 * t249 - t109 * t188 + t139 * t149 + t247 * t82;
t38 = t159 + t281 * t222 + (-t182 - t280) * t221;
t35 = -t270 * t208 - t273 * t211 + t271 * t216 - t280 * t303;
t34 = t270 * t206 + t273 * t209 - t272 * t216 + t281 * t303;
t31 = t221 * t58 + t222 * t59 + t226 * t75;
t30 = t209 * t58 + t211 * t59 - t303 * t75;
t29 = -t280 * t209 + t281 * t211 + t298;
t28 = -t247 * t75 - t248 * t59 - t249 * t58;
t27 = t221 * t53 + t222 * t54 + t226 * t66;
t26 = t221 * t51 + t222 * t52 + t226 * t65;
t25 = t209 * t53 + t211 * t54 - t303 * t66;
t24 = t209 * t51 + t211 * t52 - t303 * t65;
t23 = -t247 * t66 - t248 * t54 - t249 * t53;
t22 = -t247 * t65 - t248 * t52 - t249 * t51;
t15 = t221 * t47 + t222 * t48 + t226 * t56;
t14 = t221 * t43 + t222 * t44 + t226 * t50;
t13 = t221 * t41 + t222 * t42 + t226 * t49;
t9 = t20 * t221 + t21 * t222 + t226 * t36;
t8 = t18 * t221 + t19 * t222 + t226 * t33;
t7 = t16 * t221 + t17 * t222 + t226 * t32;
t5 = t149 * t58 + t151 * t59 + t188 * t75 - t20 * t249 - t21 * t248 - t247 * t36;
t2 = t149 * t53 + t151 * t54 - t18 * t249 + t188 * t66 - t19 * t248 - t247 * t33;
t1 = t149 * t51 + t151 * t52 - t16 * t249 - t17 * t248 + t188 * t65 - t247 * t32;
t10 = [0; 0; 0; (m(4) * t179 + m(5) * t121 + t256 * t297) * t222 + (-m(4) * t180 - m(5) * t122 - t255 * t297) * t221 + 0.2e1 * (m(5) / 0.2e1 + t297) * (-t221 * t182 + t159); (m(4) * t125 + m(5) * t67 + m(6) * t38) * t237 + ((-m(4) * t138 - m(5) * t85 - m(6) * t46) * t236 + (m(4) * t137 + m(5) * t84 + m(6) * t45) * t234) * t235; (t38 * t55 + t45 * t61 + t46 * t62) * t259 + t222 * t8 + t221 * t7 + t226 * t9 + t221 * t13 + t226 * t15 + t222 * t14 + (t67 * t88 + t84 * t93 + t85 * t94) * t261 + t226 * ((-t161 * t216 + t163 * t215 + t174 * t226 + t176 * t303 + t178 * t220) * t222 + (-t160 * t216 + t162 * t215 + t173 * t226 + t175 * t303 + t177 * t220) * t221 + (-t196 * t216 + t197 * t215 + t199 * t226 + t200 * t303 + t201 * t220) * t226) + t221 * ((-t161 * t206 - t163 * t205 + t174 * t221 - t176 * t209 + t178 * t210) * t222 + (-t160 * t206 - t162 * t205 + t173 * t221 - t175 * t209 + t177 * t210) * t221 + (-t196 * t206 - t197 * t205 + t199 * t221 - t200 * t209 + t201 * t210) * t226) + t222 * ((-t161 * t208 - t163 * t207 + t174 * t222 - t176 * t211 + t178 * t212) * t222 + (-t160 * t208 - t162 * t207 + t173 * t222 - t175 * t211 + t177 * t212) * t221 + (-t196 * t208 - t197 * t207 + t199 * t222 - t200 * t211 + t201 * t212) * t226) + 0.2e1 * m(4) * ((-t164 * t226 + t198 * t221) * t137 + (t165 * t226 - t198 * t222) * t138 + (t164 * t222 - t165 * t221) * t125); t304 + (-t209 * t255 + t211 * t256 + 0.2e1 * t298) * t297; (m(6) * t29 + t304) * t237 + ((-m(5) * t64 - m(6) * t35) * t236 + (m(5) * t63 + m(6) * t34) * t234) * t235; (t29 * t55 + t34 * t61 + t35 * t62 + t38 * t60 + t45 * t68 + t46 * t69) * m(6) + (t101 * t84 + t102 * t85 + t57 * t88 + t63 * t93 + t64 * t94 + t67 * t95) * m(5) - (t9 / 0.2e1 + t15 / 0.2e1) * t303 + (t8 / 0.2e1 + t14 / 0.2e1) * t211 + (t7 / 0.2e1 + t13 / 0.2e1) * t209 + (t221 * t71 + t222 * t72 + t226 * t91 + t26) * t291 + (t221 * t73 + t222 * t74 + t226 * t92 + t27) * t290 + (t221 * t86 + t222 * t87 + t226 * t96 + t31) * t288 + t302 * t287 + t301 * t286 + t300 * t285; (t29 * t60 + t34 * t68 + t35 * t69) * t259 + t216 * t30 + (t101 * t63 + t102 * t64 + t57 * t95) * t261 - (t282 + t300) * t303 + (t216 * t87 + t301) * t211 + (t216 * t86 + t302) * t209 + (t209 * t73 + t211 * t74 - t303 * t92 + t25) * t208 + (t209 * t71 + t211 * t72 - t303 * t91 + t24) * t206; t37 * m(6); 0.2e1 * (t37 * t237 + (t234 * t39 - t236 * t40) * t235) * t297; (t37 * t55 + t38 * t70 + t39 * t61 + t40 * t62 + t45 * t89 + t46 * t90) * m(6) + t27 * t295 + t8 * t292 + t1 * t287 + t5 * t285 + t31 * t294 + t9 * t289 + t26 * t296 + t7 * t293 + t2 * t286; t23 * t290 + t211 * t2 / 0.2e1 + t22 * t291 + t209 * t1 / 0.2e1 + (t29 * t70 + t34 * t89 + t35 * t90 + t37 * t60 + t39 * t68 + t40 * t69) * m(6) + t30 * t294 + t6 * t289 + t28 * t288 - t303 * t5 / 0.2e1 + t25 * t295 + t4 * t292 + t24 * t296 + t3 * t293; (t37 * t70 + t39 * t89 + t90 * t40) * t259 + t151 * t23 - t248 * t2 + t149 * t22 - t249 * t1 + t188 * t28 - t247 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
