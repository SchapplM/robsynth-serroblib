% Calculate time derivative of joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:29
% EndTime: 2019-12-05 15:56:47
% DurationCPUTime: 9.08s
% Computational Cost: add. (40767->781), mult. (75783->1180), div. (0->0), fcn. (85725->12), ass. (0->305)
t268 = sin(pkin(9));
t271 = cos(pkin(9));
t277 = cos(qJ(2));
t272 = cos(pkin(5));
t275 = sin(qJ(2));
t322 = t272 * t275;
t256 = t268 * t277 + t271 * t322;
t246 = t256 * qJD(2);
t301 = t268 * t322;
t305 = qJD(2) * t277;
t248 = -qJD(2) * t301 + t271 * t305;
t321 = t272 * t277;
t257 = t268 * t321 + t271 * t275;
t269 = sin(pkin(5));
t281 = -t268 * t275 + t271 * t321;
t304 = pkin(10) + qJ(4);
t266 = sin(t304);
t298 = cos(t304);
t325 = t269 * t271;
t226 = t256 * t298 - t266 * t325;
t274 = sin(qJ(5));
t276 = cos(qJ(5));
t189 = -t226 * t274 - t276 * t281;
t190 = t226 * t276 - t274 * t281;
t289 = t269 * t298;
t279 = -t256 * t266 - t271 * t289;
t109 = Icges(6,5) * t190 + Icges(6,6) * t189 - Icges(6,3) * t279;
t111 = Icges(6,4) * t190 + Icges(6,2) * t189 - Icges(6,6) * t279;
t113 = Icges(6,1) * t190 + Icges(6,4) * t189 - Icges(6,5) * t279;
t245 = t281 * qJD(2);
t184 = qJD(4) * t279 + t245 * t298;
t125 = -qJD(5) * t190 - t184 * t274 + t246 * t276;
t126 = qJD(5) * t189 + t184 * t276 + t246 * t274;
t183 = qJD(4) * t226 + t245 * t266;
t75 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t183;
t77 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t183;
t79 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t183;
t16 = t109 * t183 + t111 * t125 + t113 * t126 + t189 * t77 + t190 * t79 - t279 * t75;
t258 = t271 * t277 - t301;
t326 = t268 * t269;
t228 = t258 * t298 + t266 * t326;
t191 = -t228 * t274 + t257 * t276;
t192 = t228 * t276 + t257 * t274;
t280 = -t258 * t266 + t268 * t289;
t110 = Icges(6,5) * t192 + Icges(6,6) * t191 - Icges(6,3) * t280;
t112 = Icges(6,4) * t192 + Icges(6,2) * t191 - Icges(6,6) * t280;
t114 = Icges(6,1) * t192 + Icges(6,4) * t191 - Icges(6,5) * t280;
t247 = t257 * qJD(2);
t186 = qJD(4) * t280 - t247 * t298;
t127 = -qJD(5) * t192 - t186 * t274 + t248 * t276;
t128 = qJD(5) * t191 + t186 * t276 + t248 * t274;
t185 = qJD(4) * t228 - t247 * t266;
t76 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t185;
t78 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t185;
t80 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t185;
t17 = t110 * t183 + t112 * t125 + t114 * t126 + t189 * t78 + t190 * t80 - t279 * t76;
t306 = qJD(2) * t275;
t324 = t269 * t275;
t242 = t266 * t324 - t272 * t298;
t224 = -qJD(4) * t242 + t289 * t305;
t243 = t266 * t272 + t275 * t289;
t323 = t269 * t277;
t282 = -t243 * t276 + t274 * t323;
t299 = t269 * t306;
t154 = qJD(5) * t282 - t224 * t274 + t276 * t299;
t232 = -t243 * t274 - t276 * t323;
t155 = qJD(5) * t232 + t224 * t276 + t274 * t299;
t223 = t266 * t269 * t305 + qJD(4) * t243;
t100 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t223;
t137 = -Icges(6,5) * t282 + Icges(6,6) * t232 + Icges(6,3) * t242;
t138 = -Icges(6,4) * t282 + Icges(6,2) * t232 + Icges(6,6) * t242;
t139 = -Icges(6,1) * t282 + Icges(6,4) * t232 + Icges(6,5) * t242;
t98 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t223;
t99 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t223;
t32 = t100 * t190 + t125 * t138 + t126 * t139 + t137 * t183 + t189 * t99 - t279 * t98;
t53 = -t109 * t279 + t111 * t189 + t113 * t190;
t54 = -t110 * t279 + t112 * t189 + t114 * t190;
t67 = -t137 * t279 + t138 * t189 + t139 * t190;
t3 = -t16 * t281 + t17 * t257 + t53 * t246 + t54 * t248 + (-t277 * t32 + t306 * t67) * t269;
t117 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t246;
t119 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t246;
t121 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t246;
t142 = Icges(5,5) * t226 + Icges(5,6) * t279 - Icges(5,3) * t281;
t144 = Icges(5,4) * t226 + Icges(5,2) * t279 - Icges(5,6) * t281;
t146 = Icges(5,1) * t226 + Icges(5,4) * t279 - Icges(5,5) * t281;
t43 = -t117 * t281 + t119 * t279 + t121 * t226 + t142 * t246 - t144 * t183 + t146 * t184;
t118 = Icges(5,5) * t186 - Icges(5,6) * t185 + Icges(5,3) * t248;
t120 = Icges(5,4) * t186 - Icges(5,2) * t185 + Icges(5,6) * t248;
t122 = Icges(5,1) * t186 - Icges(5,4) * t185 + Icges(5,5) * t248;
t143 = Icges(5,5) * t228 + Icges(5,6) * t280 + Icges(5,3) * t257;
t145 = Icges(5,4) * t228 + Icges(5,2) * t280 + Icges(5,6) * t257;
t147 = Icges(5,1) * t228 + Icges(5,4) * t280 + Icges(5,5) * t257;
t44 = -t118 * t281 + t120 * t279 + t122 * t226 + t143 * t246 - t145 * t183 + t147 * t184;
t156 = Icges(5,5) * t224 - Icges(5,6) * t223 + Icges(5,3) * t299;
t157 = Icges(5,4) * t224 - Icges(5,2) * t223 + Icges(5,6) * t299;
t158 = Icges(5,1) * t224 - Icges(5,4) * t223 + Icges(5,5) * t299;
t193 = Icges(5,5) * t243 - Icges(5,6) * t242 - Icges(5,3) * t323;
t194 = Icges(5,4) * t243 - Icges(5,2) * t242 - Icges(5,6) * t323;
t195 = Icges(5,1) * t243 - Icges(5,4) * t242 - Icges(5,5) * t323;
t50 = -t156 * t281 + t157 * t279 + t158 * t226 - t183 * t194 + t184 * t195 + t193 * t246;
t85 = -t142 * t281 + t144 * t279 + t146 * t226;
t86 = -t143 * t281 + t145 * t279 + t147 * t226;
t95 = -t193 * t281 + t194 * t279 + t195 * t226;
t347 = t85 * t246 + t86 * t248 - t43 * t281 + t44 * t257 + (-t277 * t50 + t306 * t95) * t269 + t3;
t18 = t109 * t185 + t111 * t127 + t113 * t128 + t191 * t77 + t192 * t79 - t280 * t75;
t19 = t110 * t185 + t112 * t127 + t114 * t128 + t191 * t78 + t192 * t80 - t280 * t76;
t33 = t100 * t192 + t127 * t138 + t128 * t139 + t137 * t185 + t191 * t99 - t280 * t98;
t55 = -t109 * t280 + t111 * t191 + t113 * t192;
t56 = -t110 * t280 + t112 * t191 + t114 * t192;
t68 = -t137 * t280 + t138 * t191 + t139 * t192;
t4 = -t18 * t281 + t19 * t257 + t55 * t246 + t56 * t248 + (-t277 * t33 + t306 * t68) * t269;
t45 = t117 * t257 + t119 * t280 + t121 * t228 + t142 * t248 - t144 * t185 + t146 * t186;
t46 = t118 * t257 + t120 * t280 + t122 * t228 + t143 * t248 - t145 * t185 + t147 * t186;
t51 = t156 * t257 + t157 * t280 + t158 * t228 - t185 * t194 + t186 * t195 + t193 * t248;
t87 = t142 * t257 + t144 * t280 + t146 * t228;
t88 = t143 * t257 + t145 * t280 + t147 * t228;
t96 = t193 * t257 + t194 * t280 + t195 * t228;
t346 = t87 * t246 + t88 * t248 - t45 * t281 + t46 * t257 + (-t277 * t51 + t306 * t96) * t269 + t4;
t345 = 2 * m(5);
t344 = 2 * m(6);
t343 = t183 / 0.2e1;
t342 = t185 / 0.2e1;
t341 = t223 / 0.2e1;
t340 = -t279 / 0.2e1;
t339 = -t280 / 0.2e1;
t338 = t242 / 0.2e1;
t337 = t246 / 0.2e1;
t336 = t248 / 0.2e1;
t335 = t268 / 0.2e1;
t334 = -t271 / 0.2e1;
t333 = t272 / 0.2e1;
t270 = cos(pkin(10));
t332 = pkin(3) * t270;
t81 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t183;
t331 = pkin(4) * t184 + pkin(8) * t183 + t81;
t82 = rSges(6,1) * t128 + rSges(6,2) * t127 + rSges(6,3) * t185;
t330 = pkin(4) * t186 + pkin(8) * t185 + t82;
t329 = Icges(3,4) * t275;
t328 = Icges(3,4) * t277;
t267 = sin(pkin(10));
t327 = t267 * t272;
t101 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t223;
t319 = pkin(4) * t224 + pkin(8) * t223 + t101;
t115 = rSges(6,1) * t190 + rSges(6,2) * t189 - rSges(6,3) * t279;
t318 = pkin(4) * t226 - pkin(8) * t279 + t115;
t116 = rSges(6,1) * t192 + rSges(6,2) * t191 - rSges(6,3) * t280;
t317 = pkin(4) * t228 - pkin(8) * t280 + t116;
t149 = pkin(7) * t248 - t247 * t332;
t188 = -pkin(2) * t247 + qJ(3) * t248 + qJD(3) * t257;
t182 = t272 * t188;
t316 = t149 * t272 + t182;
t303 = t267 * t326;
t153 = pkin(3) * t303 + pkin(7) * t257 + t258 * t332;
t222 = pkin(2) * t258 + qJ(3) * t257;
t220 = t272 * t222;
t315 = t153 * t272 + t220;
t141 = -rSges(6,1) * t282 + rSges(6,2) * t232 + rSges(6,3) * t242;
t314 = pkin(4) * t243 + pkin(8) * t242 + t141;
t148 = pkin(7) * t246 + t245 * t332;
t187 = pkin(2) * t245 + qJ(3) * t246 - qJD(3) * t281;
t313 = -t148 - t187;
t302 = t267 * t325;
t152 = -pkin(3) * t302 - pkin(7) * t281 + t256 * t332;
t221 = pkin(2) * t256 - qJ(3) * t281;
t312 = -t152 - t221;
t311 = t187 * t326 + t188 * t325;
t259 = (pkin(2) * t275 - qJ(3) * t277) * t269;
t310 = -pkin(3) * t327 - (-pkin(7) * t277 + t275 * t332) * t269 - t259;
t309 = t221 * t326 + t222 * t325;
t239 = (-qJD(3) * t277 + (pkin(2) * t277 + qJ(3) * t275) * qJD(2)) * t269;
t307 = qJD(2) * t269;
t308 = -(pkin(7) * t275 + t277 * t332) * t307 - t239;
t253 = -t267 * t324 + t270 * t272;
t254 = t270 * t324 + t327;
t297 = (-rSges(4,1) * t254 - rSges(4,2) * t253 + rSges(4,3) * t323 - t259) * t269;
t291 = rSges(4,1) * t270 - rSges(4,2) * t267;
t296 = (-(rSges(4,3) * t275 + t277 * t291) * t307 - t239) * t269;
t295 = t148 * t326 + t149 * t325 + t311;
t294 = t152 * t326 + t153 * t325 + t309;
t159 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t299;
t293 = (-t159 + t308) * t269;
t196 = rSges(5,1) * t243 - rSges(5,2) * t242 - rSges(5,3) * t323;
t292 = (-t196 + t310) * t269;
t288 = Icges(4,1) * t270 - Icges(4,4) * t267;
t287 = Icges(4,4) * t270 - Icges(4,2) * t267;
t286 = Icges(4,5) * t270 - Icges(4,6) * t267;
t285 = -(Icges(4,4) * t254 + Icges(4,2) * t253 - Icges(4,6) * t323) * t267 + (Icges(4,1) * t254 + Icges(4,4) * t253 - Icges(4,5) * t323) * t270;
t284 = (t308 - t319) * t269;
t283 = (t310 - t314) * t269;
t235 = -t256 * t267 - t270 * t325;
t236 = t256 * t270 - t302;
t237 = -t258 * t267 + t270 * t326;
t238 = t258 * t270 + t303;
t278 = (-(Icges(4,4) * t238 + Icges(4,2) * t237 + Icges(4,6) * t257) * t267 + (Icges(4,1) * t238 + Icges(4,4) * t237 + Icges(4,5) * t257) * t270) * t268 - (-(Icges(4,4) * t236 + Icges(4,2) * t235 - Icges(4,6) * t281) * t267 + (Icges(4,1) * t236 + Icges(4,4) * t235 - Icges(4,5) * t281) * t270) * t271;
t252 = (rSges(3,1) * t277 - rSges(3,2) * t275) * t307;
t251 = (Icges(3,1) * t277 - t329) * t307;
t250 = (-Icges(3,2) * t275 + t328) * t307;
t249 = (Icges(3,5) * t277 - Icges(3,6) * t275) * t307;
t244 = t272 * rSges(3,3) + (rSges(3,1) * t275 + rSges(3,2) * t277) * t269;
t241 = Icges(3,5) * t272 + (Icges(3,1) * t275 + t328) * t269;
t240 = Icges(3,6) * t272 + (Icges(3,2) * t277 + t329) * t269;
t231 = (Icges(4,5) * t275 + t277 * t288) * t307;
t230 = (Icges(4,6) * t275 + t277 * t287) * t307;
t229 = (Icges(4,3) * t275 + t277 * t286) * t307;
t219 = -rSges(3,1) * t247 - rSges(3,2) * t248;
t218 = rSges(3,1) * t245 - rSges(3,2) * t246;
t216 = -Icges(3,1) * t247 - Icges(3,4) * t248;
t215 = Icges(3,1) * t245 - Icges(3,4) * t246;
t214 = -Icges(3,4) * t247 - Icges(3,2) * t248;
t213 = Icges(3,4) * t245 - Icges(3,2) * t246;
t212 = -Icges(3,5) * t247 - Icges(3,6) * t248;
t211 = Icges(3,5) * t245 - Icges(3,6) * t246;
t208 = rSges(3,1) * t258 - rSges(3,2) * t257 + rSges(3,3) * t326;
t207 = rSges(3,1) * t256 + rSges(3,2) * t281 - rSges(3,3) * t325;
t203 = Icges(3,1) * t258 - Icges(3,4) * t257 + Icges(3,5) * t326;
t202 = Icges(3,1) * t256 + Icges(3,4) * t281 - Icges(3,5) * t325;
t201 = Icges(3,4) * t258 - Icges(3,2) * t257 + Icges(3,6) * t326;
t200 = Icges(3,4) * t256 + Icges(3,2) * t281 - Icges(3,6) * t325;
t197 = Icges(4,5) * t254 + Icges(4,6) * t253 - Icges(4,3) * t323;
t179 = rSges(4,3) * t248 - t247 * t291;
t178 = rSges(4,3) * t246 + t245 * t291;
t177 = Icges(4,5) * t248 - t247 * t288;
t176 = Icges(4,5) * t246 + t245 * t288;
t175 = Icges(4,6) * t248 - t247 * t287;
t174 = Icges(4,6) * t246 + t245 * t287;
t173 = Icges(4,3) * t248 - t247 * t286;
t172 = Icges(4,3) * t246 + t245 * t286;
t168 = rSges(4,1) * t238 + rSges(4,2) * t237 + rSges(4,3) * t257;
t167 = rSges(4,1) * t236 + rSges(4,2) * t235 - rSges(4,3) * t281;
t161 = Icges(4,5) * t238 + Icges(4,6) * t237 + Icges(4,3) * t257;
t160 = Icges(4,5) * t236 + Icges(4,6) * t235 - Icges(4,3) * t281;
t151 = rSges(5,1) * t228 + rSges(5,2) * t280 + rSges(5,3) * t257;
t150 = rSges(5,1) * t226 + rSges(5,2) * t279 - rSges(5,3) * t281;
t131 = (t218 * t268 + t219 * t271) * t269;
t124 = rSges(5,1) * t186 - rSges(5,2) * t185 + rSges(5,3) * t248;
t123 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t246;
t108 = -t151 * t323 - t196 * t257;
t107 = t150 * t323 - t196 * t281;
t106 = (-t167 - t221) * t272 + t271 * t297;
t105 = t168 * t272 + t268 * t297 + t220;
t104 = -t193 * t323 - t194 * t242 + t195 * t243;
t103 = (-t178 - t187) * t272 + t271 * t296;
t102 = t179 * t272 + t268 * t296 + t182;
t97 = t150 * t257 + t151 * t281;
t94 = (t167 * t268 + t168 * t271) * t269 + t309;
t93 = (t178 * t268 + t179 * t271) * t269 + t311;
t92 = t116 * t242 + t141 * t280;
t91 = -t115 * t242 - t141 * t279;
t90 = -t143 * t323 - t145 * t242 + t147 * t243;
t89 = -t142 * t323 - t144 * t242 + t146 * t243;
t84 = (-t150 + t312) * t272 + t271 * t292;
t83 = t151 * t272 + t268 * t292 + t315;
t74 = t137 * t242 + t138 * t232 - t139 * t282;
t73 = -t115 * t280 + t116 * t279;
t72 = -t257 * t314 - t317 * t323;
t71 = -t281 * t314 + t318 * t323;
t70 = -t257 * t159 - t248 * t196 + (-t124 * t277 + t151 * t306) * t269;
t69 = -t281 * t159 + t246 * t196 + (t123 * t277 - t150 * t306) * t269;
t66 = (t150 * t268 + t151 * t271) * t269 + t294;
t65 = (-t123 + t313) * t272 + t271 * t293;
t64 = t124 * t272 + t268 * t293 + t316;
t63 = t257 * t318 + t281 * t317;
t62 = -t242 * t157 + t243 * t158 - t223 * t194 + t224 * t195 + (-t156 * t277 + t193 * t306) * t269;
t61 = t123 * t257 + t124 * t281 + t150 * t248 - t151 * t246;
t60 = t110 * t242 + t112 * t232 - t114 * t282;
t59 = t109 * t242 + t111 * t232 - t113 * t282;
t58 = (t312 - t318) * t272 + t271 * t283;
t57 = t268 * t283 + t272 * t317 + t315;
t52 = (t123 * t268 + t124 * t271) * t269 + t295;
t49 = (t268 * t318 + t271 * t317) * t269 + t294;
t48 = -t242 * t120 + t243 * t122 - t223 * t145 + t224 * t147 + (-t118 * t277 + t143 * t306) * t269;
t47 = -t242 * t119 + t243 * t121 - t223 * t144 + t224 * t146 + (-t117 * t277 + t142 * t306) * t269;
t42 = t101 * t280 + t116 * t223 - t141 * t185 + t242 * t82;
t41 = -t101 * t279 - t115 * t223 + t141 * t183 - t242 * t81;
t40 = (t313 - t331) * t272 + t271 * t284;
t39 = t268 * t284 + t272 * t330 + t316;
t38 = -t100 * t282 + t137 * t223 + t138 * t154 + t139 * t155 + t232 * t99 + t242 * t98;
t37 = t115 * t185 - t116 * t183 + t279 * t82 - t280 * t81;
t36 = -t319 * t257 - t314 * t248 + (-t277 * t330 + t306 * t317) * t269;
t35 = -t319 * t281 + t314 * t246 + (t277 * t331 - t306 * t318) * t269;
t34 = (t268 * t331 + t271 * t330) * t269 + t295;
t31 = t272 * t74 + (t268 * t60 - t271 * t59) * t269;
t30 = t257 * t60 - t281 * t59 - t323 * t74;
t29 = t242 * t74 - t279 * t59 - t280 * t60;
t28 = -t246 * t317 + t248 * t318 + t257 * t331 + t281 * t330;
t27 = t272 * t68 + (t268 * t56 - t271 * t55) * t269;
t26 = t272 * t67 + (t268 * t54 - t271 * t53) * t269;
t25 = t257 * t56 - t281 * t55 - t323 * t68;
t24 = t257 * t54 - t281 * t53 - t323 * t67;
t23 = t242 * t68 - t279 * t55 - t280 * t56;
t22 = t242 * t67 - t279 * t53 - t280 * t54;
t21 = t110 * t223 + t112 * t154 + t114 * t155 + t232 * t78 + t242 * t76 - t282 * t80;
t20 = t109 * t223 + t111 * t154 + t113 * t155 + t232 * t77 + t242 * t75 - t282 * t79;
t15 = t272 * t62 + (t268 * t48 - t271 * t47) * t269;
t14 = t272 * t51 + (t268 * t46 - t271 * t45) * t269;
t13 = t272 * t50 + (t268 * t44 - t271 * t43) * t269;
t12 = t89 * t246 + t90 * t248 - t47 * t281 + t48 * t257 + (t104 * t306 - t277 * t62) * t269;
t9 = t272 * t38 + (-t20 * t271 + t21 * t268) * t269;
t8 = t272 * t33 + (-t18 * t271 + t19 * t268) * t269;
t7 = t272 * t32 + (-t16 * t271 + t17 * t268) * t269;
t6 = -t20 * t281 + t21 * t257 + t59 * t246 + t60 * t248 + (-t277 * t38 + t306 * t74) * t269;
t5 = t183 * t59 + t185 * t60 - t20 * t279 - t21 * t280 + t223 * t74 + t242 * t38;
t2 = -t18 * t279 + t183 * t55 + t185 * t56 - t19 * t280 + t223 * t68 + t242 * t33;
t1 = -t16 * t279 - t17 * t280 + t183 * t53 + t185 * t54 + t223 * t67 + t242 * t32;
t10 = [0; m(3) * t131 + m(4) * t93 + m(5) * t52 + m(6) * t34; t8 * t326 - t7 * t325 - t13 * t325 + t14 * t326 - ((-t201 * t246 + t203 * t245 - t212 * t325 + t214 * t281 + t216 * t256) * t326 - (-t200 * t246 + t202 * t245 - t211 * t325 + t213 * t281 + t215 * t256) * t325 + (-t240 * t246 + t241 * t245 - t249 * t325 + t250 * t281 + t251 * t256) * t272) * t325 - ((t246 * t197 - t229 * t281 + t235 * t230 + t236 * t231 + t245 * t285) * t272 + ((t161 * t246 - t173 * t281 + t175 * t235 + t177 * t236) * t268 - (t160 * t246 - t172 * t281 + t174 * t235 + t176 * t236) * t271 + t278 * t245) * t269) * t325 + ((-t201 * t248 - t203 * t247 + t212 * t326 - t214 * t257 + t216 * t258) * t326 - (-t200 * t248 - t202 * t247 + t211 * t326 - t213 * t257 + t215 * t258) * t325 + (-t240 * t248 - t241 * t247 + t249 * t326 - t250 * t257 + t251 * t258) * t272) * t326 + ((t248 * t197 + t257 * t229 + t237 * t230 + t238 * t231 - t247 * t285) * t272 + ((t161 * t248 + t173 * t257 + t175 * t237 + t177 * t238) * t268 - (t160 * t248 + t172 * t257 + t174 * t237 + t176 * t238) * t271 - t278 * t247) * t269) * t326 + t272 * t9 + (t34 * t49 + t39 * t57 + t40 * t58) * t344 + t272 * t15 + (t52 * t66 + t64 * t83 + t65 * t84) * t345 + 0.2e1 * m(4) * (t102 * t105 + t103 * t106 + t93 * t94) + t272 * (t272 ^ 2 * t249 + (((t214 * t277 + t216 * t275) * t268 - (t213 * t277 + t215 * t275) * t271 + ((-t201 * t275 + t203 * t277) * t268 - (-t200 * t275 + t202 * t277) * t271) * qJD(2)) * t269 + (-t211 * t271 + t212 * t268 + t250 * t277 + t251 * t275 + (-t240 * t275 + t241 * t277) * qJD(2)) * t272) * t269) + t272 * ((t253 * t230 + t254 * t231) * t272 + ((t253 * t175 + t254 * t177) * t268 - (t253 * t174 + t254 * t176) * t271 + (-t229 * t272 + (t172 * t271 - t173 * t268) * t269) * t277 + ((t197 * t272 + (-t160 * t271 + t161 * t268) * t269) * t275 + (t269 * t278 + t272 * t285) * t277) * qJD(2)) * t269) + 0.2e1 * m(3) * ((-t207 * t272 - t244 * t325) * (-t218 * t272 - t252 * t325) + (t208 * t272 - t244 * t326) * (t219 * t272 - t252 * t326) + (t207 * t268 + t208 * t271) * t269 * t131); (m(4) + m(5) + m(6)) * t299; m(6) * (t246 * t57 + t248 * t58 - t281 * t39 + t257 * t40 + (-t277 * t34 + t306 * t49) * t269) + m(5) * (t246 * t83 + t248 * t84 - t281 * t64 + t257 * t65 + (-t277 * t52 + t306 * t66) * t269) + m(4) * (-t281 * t102 + t257 * t103 + t246 * t105 + t248 * t106 + (-t277 * t93 + t306 * t94) * t269); 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (-t269 ^ 2 * t275 * t305 - t246 * t281 + t248 * t257); m(5) * t61 + m(6) * t28; t27 * t336 + t26 * t337 + (t8 / 0.2e1 + t14 / 0.2e1) * t257 - (t7 / 0.2e1 + t13 / 0.2e1) * t281 + m(6) * (t28 * t49 + t34 * t63 + t35 * t58 + t36 * t57 + t39 * t72 + t40 * t71) + m(5) * (t107 * t65 + t108 * t64 + t52 * t97 + t61 * t66 + t69 * t84 + t70 * t83) + (t6 / 0.2e1 + t12 / 0.2e1 + t95 * t337 + t96 * t336) * t272 + ((t268 * t86 - t271 * t85) * t337 + (t268 * t88 - t271 * t87) * t336 + (-t9 / 0.2e1 - t15 / 0.2e1) * t277 + (t31 / 0.2e1 + t104 * t333 + (t268 * t90 - t271 * t89) * t269 / 0.2e1) * t306 + t346 * t335 + t347 * t334) * t269; m(5) * (t107 * t248 + t108 * t246 - t70 * t281 + t69 * t257 + (-t277 * t61 + t306 * t97) * t269) + m(6) * (t72 * t246 + t71 * t248 - t36 * t281 + t35 * t257 + (-t277 * t28 + t306 * t63) * t269); (-t12 - t6) * t323 + t346 * t257 - t347 * t281 + (t257 * t88 - t281 * t87 - t323 * t96 + t25) * t248 + (t257 * t86 - t281 * t85 - t323 * t95 + t24) * t246 + (-t104 * t323 + t257 * t90 - t281 * t89 + t30) * t299 + (t28 * t63 + t35 * t71 + t36 * t72) * t344 + (t107 * t69 + t108 * t70 + t61 * t97) * t345; m(6) * t37; t5 * t333 + m(6) * (t34 * t73 + t37 * t49 + t39 * t92 + t40 * t91 + t41 * t58 + t42 * t57) + t27 * t342 + t8 * t339 + t31 * t341 + t9 * t338 + t26 * t343 + t7 * t340 + (t1 * t334 + t2 * t335) * t269; m(6) * (t92 * t246 + t91 * t248 - t42 * t281 + t41 * t257 + (-t277 * t37 + t306 * t73) * t269); m(6) * (t28 * t73 + t35 * t91 + t36 * t92 + t37 * t63 + t41 * t71 + t42 * t72) + t25 * t342 + t4 * t339 + t30 * t341 + t6 * t338 + t24 * t343 + t3 * t340 + t23 * t336 + t257 * t2 / 0.2e1 + t22 * t337 - t281 * t1 / 0.2e1 + (t29 * t306 / 0.2e1 - t277 * t5 / 0.2e1) * t269; (t37 * t73 + t41 * t91 + t42 * t92) * t344 + t185 * t23 - t280 * t2 + t183 * t22 - t279 * t1 + t223 * t29 + t242 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
