% Calculate time derivative of joint inertia matrix for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:35
% DurationCPUTime: 4.73s
% Computational Cost: add. (20582->586), mult. (39637->886), div. (0->0), fcn. (42680->10), ass. (0->270)
t237 = sin(qJ(3));
t239 = cos(qJ(3));
t232 = sin(pkin(9));
t259 = qJD(3) * t232;
t213 = (-Icges(4,5) * t237 - Icges(4,6) * t239) * t259;
t234 = cos(pkin(9));
t300 = t213 * t234;
t299 = 2 * m(5);
t298 = 2 * m(6);
t235 = cos(pkin(8));
t233 = sin(pkin(8));
t276 = t233 * t239;
t253 = t234 * t276;
t258 = qJD(3) * t237;
t210 = qJD(3) * t253 - t235 * t258;
t297 = t210 / 0.2e1;
t272 = t235 * t239;
t277 = t233 * t237;
t221 = t234 * t272 + t277;
t212 = t221 * qJD(3);
t296 = t212 / 0.2e1;
t218 = t234 * t277 + t272;
t295 = t218 / 0.2e1;
t273 = t235 * t237;
t220 = t234 * t273 - t276;
t294 = t220 / 0.2e1;
t293 = -t234 / 0.2e1;
t238 = cos(qJ(4));
t292 = pkin(4) * t238;
t290 = pkin(4) * qJD(4);
t219 = t253 - t273;
t231 = qJ(4) + qJ(5);
t228 = sin(t231);
t209 = t218 * qJD(3);
t230 = qJD(4) + qJD(5);
t283 = t232 * t233;
t248 = t230 * t283 - t209;
t229 = cos(t231);
t284 = t229 * t230;
t150 = -t219 * t284 - t248 * t228;
t285 = t228 * t230;
t151 = -t219 * t285 + t248 * t229;
t100 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t210;
t192 = -t219 * t228 + t229 * t283;
t193 = t219 * t229 + t228 * t283;
t130 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t218;
t289 = t220 * t100 + t212 * t130;
t211 = t220 * qJD(3);
t282 = t232 * t235;
t247 = t230 * t282 - t211;
t152 = -t221 * t284 - t247 * t228;
t153 = -t221 * t285 + t247 * t229;
t101 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t212;
t194 = -t221 * t228 + t229 * t282;
t195 = t221 * t229 + t228 * t282;
t131 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t220;
t257 = qJD(3) * t239;
t249 = t232 * t257;
t280 = t232 * t237;
t288 = t101 * t280 + t131 * t249;
t287 = Icges(4,4) * t237;
t286 = Icges(4,4) * t239;
t236 = sin(qJ(4));
t281 = t232 * t236;
t279 = t232 * t238;
t278 = t232 * t239;
t275 = t234 * t236;
t274 = t234 * t238;
t271 = t236 * t239;
t196 = -t219 * t236 + t233 * t279;
t242 = t196 * qJD(4);
t112 = pkin(4) * t242 + pkin(7) * t210 - t292 * t209;
t270 = t100 + t112;
t198 = -t221 * t236 + t235 * t279;
t241 = t198 * qJD(4);
t113 = pkin(4) * t241 + pkin(7) * t212 - t292 * t211;
t269 = -t101 - t113;
t254 = t235 * t281;
t199 = t221 * t238 + t254;
t157 = -t199 * qJD(4) + t211 * t236;
t158 = -t211 * t238 + t241;
t111 = rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t212;
t186 = -t211 * pkin(3) + t212 * pkin(6);
t268 = -t111 - t186;
t250 = t232 * t258;
t244 = t230 * t234 + t250;
t256 = t230 * t278;
t188 = t244 * t228 - t229 * t256;
t189 = -t228 * t256 - t244 * t229;
t123 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t249;
t207 = -t228 * t278 - t229 * t234;
t208 = -t228 * t234 + t229 * t278;
t162 = rSges(6,1) * t208 + rSges(6,2) * t207 + rSges(6,3) * t280;
t267 = t218 * t123 + t210 * t162;
t154 = -t274 * t290 + (-t271 * t290 + (pkin(7) * t239 - t292 * t237) * qJD(3)) * t232;
t266 = t123 + t154;
t255 = t233 * t281;
t132 = pkin(4) * t255 + pkin(7) * t218 + t292 * t219;
t265 = t130 + t132;
t133 = pkin(4) * t254 + pkin(7) * t220 + t292 * t221;
t264 = -t131 - t133;
t141 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t220;
t191 = t221 * pkin(3) + t220 * pkin(6);
t263 = -t141 - t191;
t175 = -pkin(4) * t275 + (pkin(7) * t237 + t292 * t239) * t232;
t262 = t162 + t175;
t185 = -t209 * pkin(3) + t210 * pkin(6);
t217 = (-pkin(3) * t237 + pkin(6) * t239) * t259;
t261 = t234 * t185 + t217 * t283;
t190 = t219 * pkin(3) + t218 * pkin(6);
t224 = (pkin(3) * t239 + pkin(6) * t237) * t232;
t260 = t234 * t190 + t224 * t283;
t252 = -t186 + t269;
t251 = -t191 + t264;
t124 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t218;
t126 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t218;
t128 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t218;
t94 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t210;
t96 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t210;
t98 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t210;
t27 = t124 * t210 + t126 * t150 + t128 * t151 + t192 * t96 + t193 * t98 + t218 * t94;
t125 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t220;
t127 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t220;
t129 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t220;
t95 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t212;
t97 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t212;
t99 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t212;
t28 = t125 * t210 + t127 * t150 + t129 * t151 + t192 * t97 + t193 * t99 + t218 * t95;
t120 = Icges(6,5) * t189 + Icges(6,6) * t188 + Icges(6,3) * t249;
t121 = Icges(6,4) * t189 + Icges(6,2) * t188 + Icges(6,6) * t249;
t122 = Icges(6,1) * t189 + Icges(6,4) * t188 + Icges(6,5) * t249;
t159 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t280;
t160 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t280;
t161 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t280;
t44 = t120 * t218 + t121 * t192 + t122 * t193 + t150 * t160 + t151 * t161 + t159 * t210;
t63 = t124 * t218 + t126 * t192 + t128 * t193;
t64 = t125 * t218 + t127 * t192 + t129 * t193;
t80 = t159 * t218 + t160 * t192 + t161 * t193;
t5 = t210 * t63 + t212 * t64 + t218 * t27 + t220 * t28 + (t237 * t44 + t80 * t257) * t232;
t29 = t124 * t212 + t126 * t152 + t128 * t153 + t194 * t96 + t195 * t98 + t220 * t94;
t30 = t125 * t212 + t127 * t152 + t129 * t153 + t194 * t97 + t195 * t99 + t220 * t95;
t45 = t120 * t220 + t121 * t194 + t122 * t195 + t152 * t160 + t153 * t161 + t159 * t212;
t65 = t124 * t220 + t126 * t194 + t128 * t195;
t66 = t125 * t220 + t127 * t194 + t129 * t195;
t81 = t159 * t220 + t160 * t194 + t161 * t195;
t6 = t210 * t65 + t212 * t66 + t218 * t29 + t220 * t30 + (t237 * t45 + t81 * t257) * t232;
t73 = t124 * t280 + t126 * t207 + t128 * t208;
t74 = t125 * t280 + t127 * t207 + t129 * t208;
t85 = t159 * t280 + t160 * t207 + t161 * t208;
t32 = t126 * t188 + t128 * t189 + t207 * t96 + t208 * t98 + (t124 * t257 + t237 * t94) * t232;
t33 = t127 * t188 + t129 * t189 + t207 * t97 + t208 * t99 + (t125 * t257 + t237 * t95) * t232;
t46 = t121 * t207 + t122 * t208 + t160 * t188 + t161 * t189 + (t120 * t237 + t159 * t257) * t232;
t9 = t210 * t73 + t212 * t74 + t218 * t32 + t220 * t33 + (t237 * t46 + t85 * t257) * t232;
t246 = t218 * t5 + t220 * t6 + t9 * t280 + t210 * (t218 * t63 + t220 * t64 + t80 * t280) + t212 * (t218 * t65 + t220 * t66 + t81 * t280) + (t218 * t73 + t220 * t74 + t85 * t280) * t249;
t183 = -rSges(4,1) * t209 - rSges(4,2) * t210;
t216 = (-rSges(4,1) * t237 - rSges(4,2) * t239) * t259;
t145 = t183 * t234 + t216 * t283;
t184 = -rSges(4,1) * t211 - rSges(4,2) * t212;
t146 = -t184 * t234 - t216 * t282;
t245 = t145 * t233 - t146 * t235;
t197 = t219 * t238 + t255;
t223 = t238 * t278 - t275;
t222 = -t232 * t271 - t274;
t15 = -t234 * t44 + (t233 * t27 + t235 * t28) * t232;
t16 = -t234 * t45 + (t233 * t29 + t235 * t30) * t232;
t18 = -t234 * t46 + (t233 * t32 + t235 * t33) * t232;
t243 = t15 * t295 + t16 * t294 + t18 * t280 / 0.2e1 + (-t234 * t80 + (t233 * t63 + t235 * t64) * t232) * t297 + (-t234 * t81 + (t233 * t65 + t235 * t66) * t232) * t296 + t5 * t283 / 0.2e1 + (-t234 * t85 + (t233 * t73 + t235 * t74) * t232) * t249 / 0.2e1 + t6 * t282 / 0.2e1 + t9 * t293;
t215 = (-Icges(4,1) * t237 - t286) * t259;
t214 = (-Icges(4,2) * t239 - t287) * t259;
t205 = -Icges(4,5) * t234 + (Icges(4,1) * t239 - t287) * t232;
t204 = -Icges(4,6) * t234 + (-Icges(4,2) * t237 + t286) * t232;
t201 = t222 * qJD(4) - t238 * t250;
t200 = -t223 * qJD(4) + t236 * t250;
t182 = -Icges(4,1) * t211 - Icges(4,4) * t212;
t181 = -Icges(4,1) * t209 - Icges(4,4) * t210;
t180 = -Icges(4,4) * t211 - Icges(4,2) * t212;
t179 = -Icges(4,4) * t209 - Icges(4,2) * t210;
t178 = -Icges(4,5) * t211 - Icges(4,6) * t212;
t177 = -Icges(4,5) * t209 - Icges(4,6) * t210;
t176 = t190 * t282;
t174 = rSges(5,1) * t223 + rSges(5,2) * t222 + rSges(5,3) * t280;
t173 = Icges(5,1) * t223 + Icges(5,4) * t222 + Icges(5,5) * t280;
t172 = Icges(5,4) * t223 + Icges(5,2) * t222 + Icges(5,6) * t280;
t171 = Icges(5,5) * t223 + Icges(5,6) * t222 + Icges(5,3) * t280;
t170 = rSges(4,1) * t221 - rSges(4,2) * t220 + rSges(4,3) * t282;
t169 = rSges(4,1) * t219 - rSges(4,2) * t218 + rSges(4,3) * t283;
t168 = Icges(4,1) * t221 - Icges(4,4) * t220 + Icges(4,5) * t282;
t167 = Icges(4,1) * t219 - Icges(4,4) * t218 + Icges(4,5) * t283;
t166 = Icges(4,4) * t221 - Icges(4,2) * t220 + Icges(4,6) * t282;
t165 = Icges(4,4) * t219 - Icges(4,2) * t218 + Icges(4,6) * t283;
t163 = t185 * t282;
t156 = -t209 * t238 + t242;
t155 = -t197 * qJD(4) + t209 * t236;
t149 = t218 * t162;
t147 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t249;
t144 = Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t249;
t143 = Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t249;
t142 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t249;
t140 = rSges(5,1) * t197 + rSges(5,2) * t196 + rSges(5,3) * t218;
t139 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t220;
t138 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t218;
t137 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t220;
t136 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t218;
t135 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t220;
t134 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t218;
t119 = t131 * t280;
t117 = (t183 * t235 - t184 * t233) * t232;
t116 = t220 * t130;
t110 = rSges(5,1) * t156 + rSges(5,2) * t155 + rSges(5,3) * t210;
t109 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t212;
t108 = Icges(5,1) * t156 + Icges(5,4) * t155 + Icges(5,5) * t210;
t107 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t212;
t106 = Icges(5,4) * t156 + Icges(5,2) * t155 + Icges(5,6) * t210;
t105 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t212;
t104 = Icges(5,5) * t156 + Icges(5,6) * t155 + Icges(5,3) * t210;
t103 = t141 * t280 - t174 * t220;
t102 = -t140 * t280 + t174 * t218;
t92 = -t162 * t220 + t119;
t91 = -t130 * t280 + t149;
t89 = t171 * t280 + t172 * t222 + t173 * t223;
t88 = t140 * t220 - t141 * t218;
t87 = t263 * t234 + (-t174 - t224) * t282;
t86 = t140 * t234 + t174 * t283 + t260;
t84 = t171 * t220 + t172 * t198 + t173 * t199;
t83 = t171 * t218 + t172 * t196 + t173 * t197;
t82 = -t131 * t218 + t116;
t79 = t176 + (t140 * t235 + t263 * t233) * t232;
t78 = t135 * t280 + t137 * t222 + t139 * t223;
t77 = t134 * t280 + t136 * t222 + t138 * t223;
t76 = t268 * t234 + (-t147 - t217) * t282;
t75 = t110 * t234 + t147 * t283 + t261;
t72 = t133 * t280 - t262 * t220 + t119;
t71 = t175 * t218 - t265 * t280 + t149;
t70 = t135 * t220 + t137 * t198 + t139 * t199;
t69 = t134 * t220 + t136 * t198 + t138 * t199;
t68 = t135 * t218 + t137 * t196 + t139 * t197;
t67 = t134 * t218 + t136 * t196 + t138 * t197;
t62 = t251 * t234 + (-t224 - t262) * t282;
t61 = t265 * t234 + t262 * t283 + t260;
t60 = t163 + (t110 * t235 + t268 * t233) * t232;
t59 = -t147 * t220 - t174 * t212 + (t111 * t237 + t141 * t257) * t232;
t58 = t147 * t218 + t174 * t210 + (-t110 * t237 - t140 * t257) * t232;
t57 = t132 * t220 + t264 * t218 + t116;
t56 = -t123 * t220 - t162 * t212 + t288;
t55 = (-t100 * t237 - t130 * t257) * t232 + t267;
t54 = t143 * t222 + t144 * t223 + t172 * t200 + t173 * t201 + (t142 * t237 + t171 * t257) * t232;
t53 = t176 + (t251 * t233 + t265 * t235) * t232;
t52 = t110 * t220 - t111 * t218 + t140 * t212 - t141 * t210;
t51 = t252 * t234 + (-t217 - t266) * t282;
t50 = t270 * t234 + t266 * t283 + t261;
t49 = t142 * t220 + t143 * t198 + t144 * t199 + t157 * t172 + t158 * t173 + t171 * t212;
t48 = t142 * t218 + t143 * t196 + t144 * t197 + t155 * t172 + t156 * t173 + t171 * t210;
t47 = -t101 * t218 - t131 * t210 + t289;
t43 = t163 + (t252 * t233 + t270 * t235) * t232;
t42 = t107 * t222 + t109 * t223 + t137 * t200 + t139 * t201 + (t105 * t237 + t135 * t257) * t232;
t41 = t106 * t222 + t108 * t223 + t136 * t200 + t138 * t201 + (t104 * t237 + t134 * t257) * t232;
t40 = (t113 * t237 + t133 * t257) * t232 - t266 * t220 - t262 * t212 + t288;
t39 = t154 * t218 + t175 * t210 + (-t270 * t237 - t265 * t257) * t232 + t267;
t38 = t105 * t220 + t107 * t198 + t109 * t199 + t135 * t212 + t137 * t157 + t139 * t158;
t37 = t104 * t220 + t106 * t198 + t108 * t199 + t134 * t212 + t136 * t157 + t138 * t158;
t36 = t105 * t218 + t107 * t196 + t109 * t197 + t135 * t210 + t137 * t155 + t139 * t156;
t35 = t104 * t218 + t106 * t196 + t108 * t197 + t134 * t210 + t136 * t155 + t138 * t156;
t22 = t112 * t220 + t132 * t212 + t264 * t210 + t269 * t218 + t289;
t21 = -t234 * t54 + (t233 * t41 + t235 * t42) * t232;
t20 = -t234 * t49 + (t233 * t37 + t235 * t38) * t232;
t19 = -t234 * t48 + (t233 * t35 + t235 * t36) * t232;
t12 = t210 * t77 + t212 * t78 + t218 * t41 + t220 * t42 + (t237 * t54 + t89 * t257) * t232;
t11 = t210 * t69 + t70 * t212 + t218 * t37 + t220 * t38 + (t237 * t49 + t84 * t257) * t232;
t10 = t67 * t210 + t212 * t68 + t218 * t35 + t220 * t36 + (t237 * t48 + t83 * t257) * t232;
t1 = [0; 0; 0; m(4) * t117 + m(5) * t60 + m(6) * t43; m(4) * t245 + m(5) * (t233 * t75 - t235 * t76) + m(6) * (t233 * t50 - t235 * t51); -t234 * t18 + (t43 * t53 + t50 * t61 + t51 * t62) * t298 + (t60 * t79 + t75 * t86 + t76 * t87) * t299 - t234 * t21 + 0.2e1 * m(4) * ((t145 * t169 - t146 * t170) * t234 + ((t169 * t235 - t170 * t233) * t117 + t245 * (-rSges(4,3) * t234 + (rSges(4,1) * t239 - rSges(4,2) * t237) * t232)) * t232) - t234 * (t234 ^ 2 * t213 + (((-t180 * t237 + t182 * t239) * t235 + (-t179 * t237 + t181 * t239) * t233 + ((-t166 * t239 - t168 * t237) * t235 + (-t165 * t239 - t167 * t237) * t233) * qJD(3)) * t232 + (-t177 * t233 - t178 * t235 + t214 * t237 - t215 * t239 + (t204 * t239 + t205 * t237) * qJD(3)) * t234) * t232) + (t15 + t19 - (-t210 * t204 - t209 * t205 - t218 * t214 + t219 * t215) * t234 + (-t165 * t210 - t167 * t209 + t177 * t283 - t179 * t218 + t181 * t219 - t300) * t283) * t283 + (t16 + t20 - (-t212 * t204 - t211 * t205 - t220 * t214 + t221 * t215) * t234 + (-t166 * t212 - t168 * t211 + t178 * t282 - t180 * t220 + t182 * t221 - t300) * t282 + (-t165 * t212 - t166 * t210 - t167 * t211 - t168 * t209 + t177 * t282 + t178 * t283 - t179 * t220 - t180 * t218 + t181 * t221 + t182 * t219) * t283) * t282; m(5) * t52 + m(6) * t22; m(5) * (t233 * t58 - t235 * t59) + m(6) * (t233 * t39 - t235 * t40); (-t12 / 0.2e1 - t210 * t83 / 0.2e1 - t212 * t84 / 0.2e1) * t234 + (t235 * t11 / 0.2e1 + t233 * t10 / 0.2e1 + (t233 * t67 + t235 * t68) * t297 + (t233 * t69 + t235 * t70) * t296 + t237 * t21 / 0.2e1 + (t89 * t293 + (t233 * t77 + t235 * t78) * t232 / 0.2e1) * t257) * t232 + m(6) * (t22 * t53 + t39 * t61 + t40 * t62 + t43 * t57 + t50 * t71 + t51 * t72) + m(5) * (t102 * t75 + t103 * t76 + t52 * t79 + t58 * t86 + t59 * t87 + t60 * t88) + t20 * t294 + t19 * t295 + t243; (t218 * t77 + t220 * t78 + t89 * t280) * t249 + t12 * t280 + t212 * (t218 * t69 + t220 * t70 + t84 * t280) + t220 * t11 + t210 * (t218 * t67 + t220 * t68 + t83 * t280) + t218 * t10 + (t22 * t57 + t39 * t71 + t40 * t72) * t298 + (t102 * t58 + t103 * t59 + t52 * t88) * t299 + t246; m(6) * t47; m(6) * (t233 * t55 - t235 * t56); m(6) * (t43 * t82 + t47 * t53 + t50 * t91 + t51 * t92 + t55 * t61 + t56 * t62) + t243; m(6) * (t22 * t82 + t39 * t91 + t40 * t92 + t47 * t57 + t55 * t71 + t56 * t72) + t246; (t47 * t82 + t55 * t91 + t56 * t92) * t298 + t246;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
