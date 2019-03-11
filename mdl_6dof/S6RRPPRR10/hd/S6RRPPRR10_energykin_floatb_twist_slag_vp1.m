% Calculate kinetic energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:00
% EndTime: 2019-03-09 09:34:05
% DurationCPUTime: 5.02s
% Computational Cost: add. (1642->365), mult. (2247->526), div. (0->0), fcn. (2141->10), ass. (0->176)
t349 = Icges(3,4) + Icges(4,6);
t348 = Icges(3,1) + Icges(4,2);
t347 = Icges(3,2) + Icges(4,3);
t268 = cos(qJ(2));
t346 = t349 * t268;
t266 = sin(qJ(2));
t345 = t349 * t266;
t344 = Icges(4,4) - Icges(3,5);
t343 = Icges(4,5) - Icges(3,6);
t342 = t347 * t266 - t346;
t341 = t348 * t268 - t345;
t340 = Icges(4,1) + Icges(3,3);
t267 = sin(qJ(1));
t269 = cos(qJ(1));
t339 = t342 * t267 - t343 * t269;
t338 = t343 * t267 + t342 * t269;
t337 = -t344 * t267 + t341 * t269;
t336 = t341 * t267 + t344 * t269;
t335 = -t347 * t268 - t345;
t334 = t348 * t266 + t346;
t333 = -t343 * t266 + t344 * t268;
t240 = -qJD(2) * t269 + V_base(5);
t241 = qJD(2) * t267 + V_base(4);
t254 = V_base(6) + qJD(1);
t332 = (t266 * t335 + t268 * t334) * t254 + (t266 * t338 + t268 * t337) * t241 + (t266 * t339 + t268 * t336) * t240;
t331 = (t344 * t266 + t343 * t268) * t254 + (-t267 * t340 + t333 * t269) * t241 + (t333 * t267 + t340 * t269) * t240;
t263 = sin(pkin(10));
t327 = pkin(4) * t263;
t264 = cos(pkin(10));
t326 = t264 * pkin(4);
t325 = Icges(2,4) * t267;
t320 = qJ(4) * t266;
t319 = t266 * t269;
t262 = pkin(10) + qJ(5);
t253 = qJ(6) + t262;
t248 = sin(t253);
t318 = t267 * t248;
t249 = cos(t253);
t317 = t267 * t249;
t251 = sin(t262);
t316 = t267 * t251;
t252 = cos(t262);
t315 = t267 * t252;
t314 = t267 * t263;
t313 = t267 * t264;
t312 = t267 * t268;
t311 = t268 * t269;
t292 = pkin(2) * t268 + qJ(3) * t266;
t207 = t292 * t267;
t238 = t267 * pkin(1) - pkin(7) * t269;
t309 = -t207 - t238;
t208 = t292 * t269;
t214 = t267 * pkin(3) + qJ(4) * t311;
t308 = -t208 - t214;
t307 = pkin(5) * t252;
t305 = qJD(3) * t266;
t304 = qJD(4) * t268;
t303 = qJD(5) * t268;
t302 = qJD(6) * t268;
t301 = V_base(5) * pkin(6) + V_base(1);
t215 = -pkin(3) * t269 + qJ(4) * t312;
t298 = -t215 + t309;
t297 = pkin(5) * t251;
t233 = pkin(2) * t266 - qJ(3) * t268;
t296 = -t233 - t320;
t206 = t269 * t303 + t241;
t232 = qJD(5) * t266 + t254;
t295 = t240 * t233 + t269 * t305 + t301;
t294 = rSges(3,1) * t268 - rSges(3,2) * t266;
t293 = -rSges(4,2) * t268 + rSges(4,3) * t266;
t205 = t267 * t303 + t240;
t239 = pkin(1) * t269 + t267 * pkin(7);
t285 = -V_base(4) * pkin(6) + t254 * t239 + V_base(2);
t284 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t283 = t240 * t320 + t269 * t304 + t295;
t282 = pkin(8) * t268 + t266 * t327;
t279 = t254 * t208 + t267 * t305 + t285;
t278 = pkin(9) * t268 + t266 * t297;
t277 = -qJD(3) * t268 + t241 * t207 + t284;
t276 = t254 * t214 + t267 * t304 + t279;
t275 = qJD(4) * t266 + t241 * t215 + t277;
t151 = t267 * t282 - t269 * t326;
t198 = pkin(8) * t266 - t268 * t327;
t274 = t240 * t198 + (-t151 + t298) * t254 + t283;
t150 = t267 * t326 + t269 * t282;
t273 = t254 * t150 + (-t198 + t296) * t241 + t276;
t272 = t241 * t151 + (-t150 + t308) * t240 + t275;
t259 = Icges(2,4) * t269;
t237 = rSges(2,1) * t269 - t267 * rSges(2,2);
t236 = t267 * rSges(2,1) + rSges(2,2) * t269;
t235 = rSges(3,1) * t266 + rSges(3,2) * t268;
t234 = -rSges(4,2) * t266 - rSges(4,3) * t268;
t231 = Icges(2,1) * t269 - t325;
t230 = Icges(2,1) * t267 + t259;
t228 = -Icges(2,2) * t267 + t259;
t227 = Icges(2,2) * t269 + t325;
t219 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t218 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t217 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t210 = qJD(6) * t266 + t232;
t204 = -t264 * t269 + t266 * t314;
t203 = t263 * t269 + t266 * t313;
t202 = t263 * t319 + t313;
t201 = t264 * t319 - t314;
t197 = -rSges(4,1) * t269 + t267 * t293;
t196 = t267 * rSges(4,1) + t269 * t293;
t195 = t267 * rSges(3,3) + t269 * t294;
t194 = -rSges(3,3) * t269 + t267 * t294;
t181 = -t252 * t269 + t266 * t316;
t180 = t251 * t269 + t266 * t315;
t179 = t251 * t319 + t315;
t178 = t252 * t319 - t316;
t177 = rSges(5,3) * t266 + (-rSges(5,1) * t263 - rSges(5,2) * t264) * t268;
t175 = Icges(5,5) * t266 + (-Icges(5,1) * t263 - Icges(5,4) * t264) * t268;
t174 = Icges(5,6) * t266 + (-Icges(5,4) * t263 - Icges(5,2) * t264) * t268;
t173 = Icges(5,3) * t266 + (-Icges(5,5) * t263 - Icges(5,6) * t264) * t268;
t171 = -t249 * t269 + t266 * t318;
t170 = t248 * t269 + t266 * t317;
t169 = t248 * t319 + t317;
t168 = t249 * t319 - t318;
t167 = t269 * t302 + t206;
t166 = t267 * t302 + t205;
t164 = rSges(6,3) * t266 + (-rSges(6,1) * t251 - rSges(6,2) * t252) * t268;
t163 = Icges(6,5) * t266 + (-Icges(6,1) * t251 - Icges(6,4) * t252) * t268;
t162 = Icges(6,6) * t266 + (-Icges(6,4) * t251 - Icges(6,2) * t252) * t268;
t161 = Icges(6,3) * t266 + (-Icges(6,5) * t251 - Icges(6,6) * t252) * t268;
t160 = rSges(7,3) * t266 + (-rSges(7,1) * t248 - rSges(7,2) * t249) * t268;
t159 = V_base(5) * rSges(2,3) - t236 * t254 + t301;
t158 = t237 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t157 = Icges(7,5) * t266 + (-Icges(7,1) * t248 - Icges(7,4) * t249) * t268;
t156 = Icges(7,6) * t266 + (-Icges(7,4) * t248 - Icges(7,2) * t249) * t268;
t155 = Icges(7,3) * t266 + (-Icges(7,5) * t248 - Icges(7,6) * t249) * t268;
t153 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t152 = pkin(9) * t266 - t268 * t297;
t149 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t312;
t148 = t202 * rSges(5,1) + t201 * rSges(5,2) + rSges(5,3) * t311;
t147 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t312;
t146 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t311;
t145 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t312;
t144 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t311;
t143 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t312;
t142 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t311;
t139 = rSges(6,1) * t181 + rSges(6,2) * t180 + rSges(6,3) * t312;
t138 = t179 * rSges(6,1) + t178 * rSges(6,2) + rSges(6,3) * t311;
t137 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t312;
t136 = Icges(6,1) * t179 + Icges(6,4) * t178 + Icges(6,5) * t311;
t135 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t312;
t134 = Icges(6,4) * t179 + Icges(6,2) * t178 + Icges(6,6) * t311;
t133 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t312;
t132 = Icges(6,5) * t179 + Icges(6,6) * t178 + Icges(6,3) * t311;
t131 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t312;
t130 = t169 * rSges(7,1) + t168 * rSges(7,2) + rSges(7,3) * t311;
t129 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t312;
t128 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t311;
t127 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t312;
t126 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t311;
t125 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t312;
t124 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t311;
t123 = t235 * t240 + (-t194 - t238) * t254 + t301;
t122 = t195 * t254 - t235 * t241 + t285;
t121 = t267 * t278 - t269 * t307;
t120 = t267 * t307 + t269 * t278;
t119 = t194 * t241 - t195 * t240 + t284;
t118 = t234 * t240 + (-t197 + t309) * t254 + t295;
t117 = t196 * t254 + (-t233 - t234) * t241 + t279;
t116 = t197 * t241 + (-t196 - t208) * t240 + t277;
t115 = t177 * t240 + (-t149 + t298) * t254 + t283;
t114 = t148 * t254 + (-t177 + t296) * t241 + t276;
t113 = t149 * t241 + (-t148 + t308) * t240 + t275;
t112 = -t139 * t232 + t164 * t205 + t274;
t111 = t138 * t232 - t164 * t206 + t273;
t110 = -t138 * t205 + t139 * t206 + t272;
t109 = -t121 * t232 - t131 * t210 + t152 * t205 + t160 * t166 + t274;
t108 = t120 * t232 + t130 * t210 - t152 * t206 - t160 * t167 + t273;
t107 = -t120 * t205 + t121 * t206 - t130 * t166 + t131 * t167 + t272;
t1 = t210 * ((t124 * t167 + t125 * t166 + t155 * t210) * t266 + ((-t126 * t249 - t128 * t248) * t167 + (-t127 * t249 - t129 * t248) * t166 + (-t156 * t249 - t157 * t248) * t210) * t268) / 0.2e1 + m(1) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + m(2) * (t153 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(3) * (t119 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t232 * ((t132 * t206 + t133 * t205 + t161 * t232) * t266 + ((-t134 * t252 - t136 * t251) * t206 + (-t135 * t252 - t137 * t251) * t205 + (-t162 * t252 - t163 * t251) * t232) * t268) / 0.2e1 + t206 * ((t132 * t311 + t178 * t134 + t179 * t136) * t206 + (t133 * t311 + t178 * t135 + t179 * t137) * t205 + (t161 * t311 + t178 * t162 + t179 * t163) * t232) / 0.2e1 + t167 * ((t124 * t311 + t168 * t126 + t169 * t128) * t167 + (t125 * t311 + t168 * t127 + t169 * t129) * t166 + (t155 * t311 + t168 * t156 + t169 * t157) * t210) / 0.2e1 + t205 * ((t132 * t312 + t134 * t180 + t136 * t181) * t206 + (t133 * t312 + t180 * t135 + t181 * t137) * t205 + (t161 * t312 + t162 * t180 + t163 * t181) * t232) / 0.2e1 + t166 * ((t124 * t312 + t126 * t170 + t128 * t171) * t167 + (t125 * t312 + t170 * t127 + t171 * t129) * t166 + (t155 * t312 + t156 * t170 + t157 * t171) * t210) / 0.2e1 + ((-t267 * t227 + t230 * t269 + Icges(1,4)) * V_base(5) + (-t267 * t228 + t269 * t231 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t269 * t227 + t267 * t230 + Icges(1,2)) * V_base(5) + (t228 * t269 + t267 * t231 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t142 * t312 + t144 * t203 + t146 * t204) * t241 + (t143 * t312 + t203 * t145 + t204 * t147) * t240 + (t173 * t312 + t174 * t203 + t175 * t204) * t254 + t331 * t269 + t332 * t267) * t240 / 0.2e1 + ((t142 * t311 + t201 * t144 + t202 * t146) * t241 + (t143 * t311 + t201 * t145 + t202 * t147) * t240 + (t173 * t311 + t201 * t174 + t202 * t175) * t254 + t332 * t269 - t331 * t267) * t241 / 0.2e1 + (((-t144 * t264 - t146 * t263 - t338) * t268 + (t142 + t337) * t266) * t241 + ((-t145 * t264 - t147 * t263 - t339) * t268 + (t143 + t336) * t266) * t240 + (Icges(2,3) + (-t174 * t264 - t175 * t263 - t335) * t268 + (t173 + t334) * t266) * t254) * t254 / 0.2e1 + t254 * V_base(4) * (Icges(2,5) * t269 - Icges(2,6) * t267) + t254 * V_base(5) * (Icges(2,5) * t267 + Icges(2,6) * t269) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
