% Calculate time derivative of joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR9_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:22
% EndTime: 2019-12-31 17:09:34
% DurationCPUTime: 7.21s
% Computational Cost: add. (9491->610), mult. (15660->890), div. (0->0), fcn. (14867->8), ass. (0->294)
t207 = sin(qJ(1));
t337 = t207 / 0.2e1;
t209 = cos(qJ(1));
t355 = -t209 / 0.2e1;
t354 = -qJD(1) / 0.2e1;
t204 = cos(pkin(7));
t191 = pkin(3) * t204 + pkin(2);
t208 = cos(qJ(2));
t306 = t208 * t209;
t205 = -pkin(6) - qJ(3);
t203 = sin(pkin(7));
t309 = t207 * t203;
t206 = sin(qJ(2));
t310 = t206 * t209;
t350 = pkin(3) * t309 - t205 * t310;
t200 = pkin(7) + qJ(4);
t194 = sin(t200);
t195 = cos(t200);
t144 = -t194 * t306 + t207 * t195;
t145 = t207 * t194 + t195 * t306;
t90 = t145 * rSges(5,1) + t144 * rSges(5,2) + rSges(5,3) * t310;
t353 = t191 * t306 + t350 + t90;
t330 = pkin(2) - t191;
t352 = t206 * t330;
t305 = qJ(3) + t205;
t351 = t208 * t305;
t322 = Icges(3,4) * t208;
t241 = -Icges(3,2) * t206 + t322;
t149 = Icges(3,6) * t207 + t209 * t241;
t323 = Icges(3,4) * t206;
t245 = Icges(3,1) * t208 - t323;
t151 = Icges(3,5) * t207 + t209 * t245;
t231 = t149 * t206 - t151 * t208;
t222 = t231 * t207;
t148 = -Icges(3,6) * t209 + t207 * t241;
t150 = -Icges(3,5) * t209 + t207 * t245;
t232 = t148 * t206 - t150 * t208;
t223 = t232 * t209;
t349 = -rSges(3,2) * t310 + t207 * rSges(3,3);
t237 = Icges(3,5) * t208 - Icges(3,6) * t206;
t146 = -Icges(3,3) * t209 + t207 * t237;
t292 = qJD(1) * t208;
t261 = -qJD(4) + t292;
t287 = qJD(2) * t209;
t267 = t206 * t287;
t348 = t207 * t261 + t267;
t289 = qJD(2) * t207;
t268 = t206 * t289;
t347 = t209 * t261 - t268;
t346 = 2 * m(3);
t345 = 2 * m(4);
t344 = 2 * m(5);
t201 = t207 ^ 2;
t202 = t209 ^ 2;
t343 = m(4) / 0.2e1;
t342 = m(5) / 0.2e1;
t239 = Icges(4,4) * t204 - Icges(4,2) * t203;
t139 = -Icges(4,6) * t208 + t206 * t239;
t341 = t139 / 0.2e1;
t243 = Icges(4,1) * t204 - Icges(4,4) * t203;
t140 = -Icges(4,5) * t208 + t206 * t243;
t340 = t140 / 0.2e1;
t339 = -t203 / 0.2e1;
t338 = t204 / 0.2e1;
t336 = -t208 / 0.2e1;
t334 = t209 / 0.2e1;
t177 = rSges(3,1) * t206 + rSges(3,2) * t208;
t333 = m(3) * t177;
t332 = pkin(2) * t208;
t197 = t207 * pkin(5);
t331 = qJD(1) / 0.2e1;
t329 = rSges(3,3) * t209;
t328 = rSges(5,3) * t206;
t327 = -rSges(4,3) - qJ(3);
t326 = -rSges(5,3) + t205;
t212 = -t206 * t305 - t208 * t330;
t314 = t203 * t209;
t284 = pkin(3) * t314;
t308 = t207 * t208;
t142 = -t194 * t308 - t195 * t209;
t143 = -t194 * t209 + t195 * t308;
t255 = -t143 * rSges(5,1) - t142 * rSges(5,2);
t311 = t206 * t207;
t89 = rSges(5,3) * t311 - t255;
t325 = t207 * t212 - t284 + t89;
t190 = pkin(2) * t306;
t165 = qJ(3) * t310 + t190;
t324 = -t165 + t353;
t321 = Icges(5,4) * t194;
t320 = Icges(5,4) * t195;
t242 = Icges(5,1) * t195 - t321;
t133 = -Icges(5,5) * t208 + t206 * t242;
t316 = t133 * t195;
t315 = t203 * t208;
t313 = t204 * t208;
t312 = t206 * t191;
t307 = t208 * t205;
t162 = -t203 * t306 + t207 * t204;
t163 = t204 * t306 + t309;
t104 = t163 * rSges(4,1) + t162 * rSges(4,2) + rSges(4,3) * t310;
t304 = -t104 - t165;
t254 = rSges(5,1) * t195 - rSges(5,2) * t194;
t135 = -rSges(5,3) * t208 + t206 * t254;
t303 = t135 + t351 - t352;
t253 = qJ(3) * t206 + t332;
t159 = qJD(2) * t253 - qJD(3) * t208;
t256 = rSges(4,1) * t204 - rSges(4,2) * t203;
t302 = -(rSges(4,3) * t206 + t208 * t256) * qJD(2) - t159;
t141 = -rSges(4,3) * t208 + t206 * t256;
t176 = pkin(2) * t206 - qJ(3) * t208;
t301 = -t141 - t176;
t164 = t253 * t207;
t300 = t207 * t164 + t209 * t165;
t293 = qJD(1) * t207;
t270 = t206 * t293;
t299 = qJD(1) * t284 + t205 * t270;
t291 = qJD(1) * t209;
t298 = rSges(3,2) * t270 + rSges(3,3) * t291;
t286 = qJD(3) * t206;
t184 = t209 * t286;
t193 = pkin(5) * t291;
t297 = t184 + t193;
t296 = t209 * pkin(1) + t197;
t295 = t201 + t202;
t147 = Icges(3,3) * t207 + t209 * t237;
t294 = qJD(1) * t147;
t290 = qJD(2) * t206;
t288 = qJD(2) * t208;
t285 = qJD(4) * t206;
t265 = t208 * t287;
t178 = qJ(3) * t265;
t183 = pkin(2) * t268;
t266 = t207 * t288;
t215 = t206 * t291 + t266;
t216 = -t207 * t292 - t267;
t283 = t164 * t291 + t207 * (qJ(3) * t215 + qJD(1) * t190 + t207 * t286 - t183) + t209 * (pkin(2) * t216 - qJ(3) * t270 + t178 + t184);
t262 = -qJD(4) * t208 + qJD(1);
t229 = t262 * t209;
t77 = t194 * t348 + t195 * t229;
t78 = t194 * t229 - t195 * t348;
t282 = t78 * rSges(5,1) + t77 * rSges(5,2) + rSges(5,3) * t265;
t161 = t204 * t308 - t314;
t225 = t203 * t308 + t204 * t209;
t97 = Icges(4,5) * t161 - Icges(4,6) * t225 + Icges(4,3) * t311;
t280 = t97 * t311;
t98 = Icges(4,5) * t163 + Icges(4,6) * t162 + Icges(4,3) * t310;
t279 = t98 * t311;
t278 = t97 * t310;
t277 = t98 * t310;
t85 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t311;
t87 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t311;
t252 = -t194 * t85 + t195 * t87;
t83 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t311;
t32 = t206 * t252 - t208 * t83;
t235 = Icges(5,5) * t195 - Icges(5,6) * t194;
t131 = -Icges(5,3) * t208 + t206 * t235;
t238 = -Icges(5,2) * t194 + t320;
t132 = -Icges(5,6) * t208 + t206 * t238;
t49 = t131 * t311 + t132 * t142 + t133 * t143;
t276 = t32 / 0.2e1 + t49 / 0.2e1;
t86 = Icges(5,4) * t145 + Icges(5,2) * t144 + Icges(5,6) * t310;
t88 = Icges(5,1) * t145 + Icges(5,4) * t144 + Icges(5,5) * t310;
t251 = -t194 * t86 + t195 * t88;
t84 = Icges(5,5) * t145 + Icges(5,6) * t144 + Icges(5,3) * t310;
t33 = t206 * t251 - t208 * t84;
t50 = t131 * t310 + t144 * t132 + t145 * t133;
t275 = t33 / 0.2e1 + t50 / 0.2e1;
t94 = (-rSges(5,1) * t194 - rSges(5,2) * t195) * t285 + (t208 * t254 + t328) * qJD(2);
t274 = -t212 * qJD(2) - t159 - t94;
t119 = qJD(1) * t225 + t203 * t267;
t120 = -qJD(1) * t161 - t204 * t267;
t272 = t120 * rSges(4,1) + t119 * rSges(4,2) + rSges(4,3) * t265;
t271 = -t176 - t303;
t269 = t132 * t288;
t264 = t288 / 0.2e1;
t263 = -t191 * t208 - pkin(1);
t114 = t301 * t209;
t71 = t271 * t209;
t230 = t262 * t207;
t79 = -t194 * t347 + t195 * t230;
t80 = t194 * t230 + t195 * t347;
t260 = t80 * rSges(5,1) + t79 * rSges(5,2);
t259 = rSges(3,1) * t208 - rSges(3,2) * t206;
t121 = qJD(1) * t162 + t203 * t268;
t122 = qJD(1) * t163 - t204 * t268;
t258 = -t122 * rSges(4,1) - t121 * rSges(4,2);
t257 = -rSges(4,1) * t161 + rSges(4,2) * t225;
t26 = t142 * t85 + t143 * t87 + t311 * t83;
t27 = t142 * t86 + t143 * t88 + t311 * t84;
t17 = t27 * t207 - t209 * t26;
t250 = t207 * t26 + t209 * t27;
t28 = t144 * t85 + t145 * t87 + t310 * t83;
t29 = t144 * t86 + t145 * t88 + t310 * t84;
t18 = t29 * t207 - t209 * t28;
t249 = t207 * t28 + t209 * t29;
t248 = t33 * t207 - t32 * t209;
t247 = t32 * t207 + t33 * t209;
t246 = -t207 * t90 + t209 * t89;
t244 = Icges(3,1) * t206 + t322;
t240 = Icges(3,2) * t208 + t323;
t236 = Icges(4,5) * t204 - Icges(4,6) * t203;
t154 = rSges(3,1) * t306 + t349;
t228 = -pkin(1) - t259;
t224 = qJD(2) * t177;
t221 = qJD(2) * t244;
t220 = qJD(2) * t240;
t219 = qJD(2) * (-Icges(3,5) * t206 - Icges(3,6) * t208);
t218 = t206 * t327 - pkin(1) - t332;
t217 = t206 * t326 + t263;
t214 = t265 - t270;
t211 = t218 * t207;
t91 = (-Icges(5,5) * t194 - Icges(5,6) * t195) * t285 + (Icges(5,3) * t206 + t208 * t235) * qJD(2);
t93 = (-Icges(5,1) * t194 - t320) * t285 + (Icges(5,5) * t206 + t208 * t242) * qJD(2);
t210 = t131 * t290 - t208 * t91 + t288 * t316 + (-t132 * t285 + t206 * t93) * t195;
t198 = t209 * pkin(5);
t170 = t259 * qJD(2);
t166 = t176 * t293;
t153 = t207 * t259 - t329;
t130 = (Icges(4,5) * t206 + t208 * t243) * qJD(2);
t129 = (Icges(4,6) * t206 + t208 * t239) * qJD(2);
t126 = t154 + t296;
t125 = t207 * t228 + t198 + t329;
t113 = t301 * t207;
t108 = t207 * t219 + t294;
t107 = -qJD(1) * t146 + t209 * t219;
t103 = rSges(4,3) * t311 - t257;
t102 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t310;
t101 = Icges(4,1) * t161 - Icges(4,4) * t225 + Icges(4,5) * t311;
t100 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t310;
t99 = Icges(4,4) * t161 - Icges(4,2) * t225 + Icges(4,6) * t311;
t92 = (-Icges(5,2) * t195 - t321) * t285 + (Icges(5,6) * t206 + t208 * t238) * qJD(2);
t82 = t177 * t289 + ((-rSges(3,3) - pkin(5)) * t207 + t228 * t209) * qJD(1);
t81 = rSges(3,1) * t216 - rSges(3,2) * t265 - pkin(1) * t293 + t193 + t298;
t73 = t296 - t304;
t72 = t198 + t211 + t257;
t70 = t271 * t207;
t69 = t207 * t147 - t209 * t231;
t68 = t207 * t146 - t223;
t67 = -t147 * t209 - t222;
t66 = -t146 * t209 - t207 * t232;
t65 = Icges(4,1) * t122 + Icges(4,4) * t121 + Icges(4,5) * t215;
t64 = Icges(4,1) * t120 + Icges(4,4) * t119 + Icges(4,5) * t214;
t63 = Icges(4,4) * t122 + Icges(4,2) * t121 + Icges(4,6) * t215;
t62 = Icges(4,4) * t120 + Icges(4,2) * t119 + Icges(4,6) * t214;
t59 = -t135 * t310 - t208 * t90;
t58 = t135 * t311 + t208 * t89;
t57 = qJD(1) * t114 + t207 * t302;
t56 = t141 * t293 + t209 * t302 + t166;
t55 = t296 + t353;
t54 = t207 * t217 + t198 + t255 + t284;
t53 = -t131 * t208 + (-t132 * t194 + t316) * t206;
t52 = t53 * t290;
t51 = t246 * t206;
t48 = t207 * t103 + t104 * t209 + t300;
t47 = t183 + (t288 * t327 - t286) * t207 + (t209 * t218 - t197) * qJD(1) + t258;
t46 = -pkin(2) * t267 + qJD(1) * t211 + t178 + t272 + t297;
t45 = rSges(5,3) * t215 + t260;
t44 = -rSges(5,3) * t270 + t282;
t43 = Icges(5,1) * t80 + Icges(5,4) * t79 + Icges(5,5) * t215;
t42 = Icges(5,1) * t78 + Icges(5,4) * t77 + Icges(5,5) * t214;
t41 = Icges(5,4) * t80 + Icges(5,2) * t79 + Icges(5,6) * t215;
t40 = Icges(5,4) * t78 + Icges(5,2) * t77 + Icges(5,6) * t214;
t39 = Icges(5,5) * t80 + Icges(5,6) * t79 + Icges(5,3) * t215;
t38 = Icges(5,5) * t78 + Icges(5,6) * t77 + Icges(5,3) * t214;
t37 = t162 * t100 + t163 * t102 + t277;
t36 = t163 * t101 + t162 * t99 + t278;
t35 = -t100 * t225 + t102 * t161 + t279;
t34 = t101 * t161 - t225 * t99 + t280;
t31 = qJD(1) * t71 + t207 * t274;
t30 = t209 * t274 + t293 * t303 + t166;
t25 = (-t286 + (t208 * t326 + t312) * qJD(2)) * t207 + ((-pkin(3) * t203 - pkin(5)) * t207 + t217 * t209) * qJD(1) - t260;
t24 = (-t307 - t312) * t287 + (t263 - t328) * t293 + t282 + t297 + t299;
t23 = t207 * t325 + t209 * t324 + t300;
t22 = (t135 * t289 + t45) * t208 + (-qJD(2) * t89 + t135 * t291 + t207 * t94) * t206;
t21 = (-t135 * t287 - t44) * t208 + (qJD(2) * t90 + t135 * t293 - t209 * t94) * t206;
t20 = (-t269 + (-qJD(4) * t133 - t92) * t206) * t194 + t210;
t19 = t207 * (rSges(4,3) * t266 - t258) + t209 * t272 + (t209 * t103 + t207 * t304) * qJD(1) + t283;
t16 = t131 * t215 + t79 * t132 + t80 * t133 + t142 * t92 + t143 * t93 + t311 * t91;
t15 = t131 * t214 + t77 * t132 + t78 * t133 + t144 * t92 + t145 * t93 + t310 * t91;
t14 = t246 * t288 + (-t207 * t44 + t209 * t45 + (-t207 * t89 - t209 * t90) * qJD(1)) * t206;
t13 = t206 * t249 - t50 * t208;
t12 = t206 * t250 - t49 * t208;
t11 = (qJD(2) * t251 - t38) * t208 + (qJD(2) * t84 - t194 * t40 + t195 * t42 + (-t194 * t88 - t195 * t86) * qJD(4)) * t206;
t10 = (qJD(2) * t252 - t39) * t208 + (qJD(2) * t83 - t194 * t41 + t195 * t43 + (-t194 * t87 - t195 * t85) * qJD(4)) * t206;
t9 = t84 * t266 + t142 * t40 + t143 * t42 + t79 * t86 + t80 * t88 + (t207 * t38 + t291 * t84) * t206;
t8 = t83 * t266 + t142 * t41 + t143 * t43 + t79 * t85 + t80 * t87 + (t207 * t39 + t291 * t83) * t206;
t7 = t84 * t265 + t144 * t40 + t145 * t42 + t77 * t86 + t78 * t88 + (t209 * t38 - t293 * t84) * t206;
t6 = t83 * t265 + t144 * t41 + t145 * t43 + t77 * t85 + t78 * t87 + (t209 * t39 - t293 * t83) * t206;
t5 = (-t178 + t44 + t299) * t209 + (t183 + t45) * t207 + (t202 * (-t307 + t352) + (-t312 - t351) * t201) * qJD(2) + (t325 * t209 + (-t165 - t324 + t350) * t207) * qJD(1) + t283;
t4 = qJD(1) * t250 + t9 * t207 - t209 * t8;
t3 = qJD(1) * t249 + t7 * t207 - t209 * t6;
t2 = (qJD(2) * t250 - t16) * t208 + (-qJD(1) * t17 + qJD(2) * t49 + t207 * t8 + t209 * t9) * t206;
t1 = (qJD(2) * t249 - t15) * t208 + (-qJD(1) * t18 + qJD(2) * t50 + t207 * t6 + t209 * t7) * t206;
t60 = [t210 + (t24 * t55 + t25 * t54) * t344 + (t46 * t73 + t47 * t72) * t345 + (t125 * t82 + t126 * t81) * t346 + (-t203 * t129 + t204 * t130) * t206 + (-Icges(4,3) * t208 + t206 * t236 - t240 + t245) * t290 + (-t133 * t285 - t206 * t92 - t269) * t194 + (-Icges(4,3) * t206 - t203 * t139 + t204 * t140 - t208 * t236 + t241 + t244) * t288; m(3) * ((-t207 * t81 - t209 * t82) * t177 + (-t125 * t209 - t126 * t207) * t170) + m(5) * (t24 * t70 + t25 * t71 + t30 * t54 + t31 * t55) + m(4) * (t113 * t46 + t114 * t47 + t56 * t72 + t57 * t73) + ((t149 * t354 + t220 * t337 + Icges(4,5) * t122 / 0.2e1 + Icges(4,6) * t121 / 0.2e1 + Icges(4,3) * t215 / 0.2e1) * t209 + (t148 * t354 + t220 * t355 - Icges(4,5) * t120 / 0.2e1 - Icges(4,6) * t119 / 0.2e1 - Icges(4,3) * t214 / 0.2e1) * t207) * t208 + ((-t126 * t333 + t162 * t341 + t163 * t340 + (-t98 / 0.2e1 + t149 / 0.2e1) * t208 + (t100 * t339 + t102 * t338 + t151 / 0.2e1) * t206 + t275) * t209 + (t125 * t333 - t225 * t341 + t161 * t340 + (-t97 / 0.2e1 + t148 / 0.2e1) * t208 + (t101 * t338 + t99 * t339 + t150 / 0.2e1) * t206 + t276) * t207) * qJD(1) + ((t201 / 0.2e1 + t202 / 0.2e1) * t237 - t222 / 0.2e1 + t223 / 0.2e1) * qJD(2) + (t11 + t15 + t119 * t139 + t120 * t140 + t162 * t129 + t163 * t130 + (-qJD(1) * t150 - t203 * t62 + t204 * t64 - t209 * t221) * t206 + (-t100 * t315 + t102 * t313 + t206 * t98) * qJD(2)) * t337 + (t10 + t16 + t121 * t139 + t122 * t140 - t225 * t129 + t161 * t130 + (qJD(1) * t151 - t203 * t63 + t204 * t65 - t207 * t221) * t206 + (t101 * t313 + t206 * t97 - t315 * t99) * qJD(2)) * t355; (t23 * t5 + t30 * t71 + t31 * t70) * t344 + t207 * t3 - t209 * t4 + (t113 * t57 + t114 * t56 + t48 * t19) * t345 - t209 * ((t209 * t108 + (t67 + t223) * qJD(1)) * t209 + (t66 * qJD(1) + (-t149 * t288 - t151 * t290 + t294) * t207 + (-t107 + (t148 * t208 + t150 * t206) * qJD(2) - t231 * qJD(1)) * t209) * t207) + t207 * ((t119 * t100 + t120 * t102 + t162 * t62 + t163 * t64 + (t36 - t279) * qJD(1)) * t207 + (-t120 * t101 - t119 * t99 - t162 * t63 - t163 * t65 + (t37 + t280) * qJD(1)) * t209) - t209 * ((-t122 * t101 - t121 * t99 + t225 * t63 - t161 * t65 + (t35 - t278) * qJD(1)) * t209 + (t121 * t100 + t122 * t102 - t225 * t62 + t161 * t64 + (t34 + t277) * qJD(1)) * t207) + t207 * ((t207 * t107 + (t68 + t222) * qJD(1)) * t207 + (t69 * qJD(1) + (t148 * t288 + t150 * t290) * t209 + (-t108 + (-t149 * t208 - t151 * t206) * qJD(2) + (t147 - t232) * qJD(1)) * t207) * t209) + ((t207 * t153 + t154 * t209) * ((qJD(1) * t153 - t209 * t224 + t298) * t209 + (-t207 * t224 + (-t154 + t349) * qJD(1)) * t207) + t295 * t177 * t170) * t346 + (t17 + (-t34 - t66) * t209 + (t35 + t67) * t207) * t293 + (t18 + (-t36 - t68) * t209 + (t37 + t69) * t207) * t291; 0.2e1 * ((t207 * t55 + t209 * t54) * t342 + (t207 * t73 + t209 * t72) * t343) * t288 + 0.2e1 * ((t207 * t24 + t209 * t25 + t291 * t55 - t293 * t54) * t342 + (t207 * t46 + t209 * t47 + t291 * t73 - t293 * t72) * t343) * t206; 0.2e1 * ((t287 * t71 + t289 * t70 - t5) * t342 + (t113 * t289 + t114 * t287 - t19) * t343) * t208 + 0.2e1 * ((qJD(2) * t23 + t207 * t31 + t209 * t30 + t291 * t70 - t293 * t71) * t342 + (qJD(2) * t48 + t113 * t291 - t114 * t293 + t207 * t57 + t209 * t56) * t343) * t206; 0.4e1 * (t343 + t342) * (-0.1e1 + t295) * t206 * t288; m(5) * (t21 * t55 + t22 * t54 + t24 * t59 + t25 * t58) + t52 + (-t20 + (t207 * t276 + t209 * t275) * qJD(2)) * t208 + ((t11 / 0.2e1 + t15 / 0.2e1) * t209 + (t10 / 0.2e1 + t16 / 0.2e1) * t207 + (-t207 * t275 + t209 * t276) * qJD(1)) * t206; m(5) * (t14 * t23 + t21 * t70 + t22 * t71 + t58 * t30 + t59 * t31 + t51 * t5) + (t18 * t264 + t13 * t331 - t2 / 0.2e1 + (qJD(1) * t33 - t10) * t336) * t209 + (t1 / 0.2e1 + t12 * t331 + t17 * t264 + (t32 * qJD(1) + t11) * t336) * t207 + (t3 * t334 + t4 * t337 + qJD(2) * t248 / 0.2e1 + (-t207 * t18 / 0.2e1 + t17 * t334) * qJD(1)) * t206; m(5) * ((-t14 + (t207 * t59 + t209 * t58) * qJD(2)) * t208 + (qJD(2) * t51 + t207 * t21 + t209 * t22 + (-t207 * t58 + t209 * t59) * qJD(1)) * t206); (t14 * t51 + t21 * t59 + t22 * t58) * t344 + (t20 * t208 - t52 + (t207 * t12 + t209 * t13 - t208 * t247) * qJD(2)) * t208 + (t209 * t1 + t207 * t2 - t208 * (t10 * t207 + t11 * t209) + (t206 * t247 - t53 * t208) * qJD(2) + (t209 * t12 - t207 * t13 + t208 * t248) * qJD(1)) * t206;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t60(1), t60(2), t60(4), t60(7); t60(2), t60(3), t60(5), t60(8); t60(4), t60(5), t60(6), t60(9); t60(7), t60(8), t60(9), t60(10);];
Mq = res;
