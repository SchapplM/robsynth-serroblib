% Calculate kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:25
% EndTime: 2019-03-09 10:03:28
% DurationCPUTime: 3.36s
% Computational Cost: add. (1144->262), mult. (2213->356), div. (0->0), fcn. (2131->6), ass. (0->133)
t327 = Icges(3,4) + Icges(4,6);
t326 = Icges(3,1) + Icges(4,2);
t325 = -Icges(3,2) - Icges(4,3);
t238 = cos(qJ(2));
t324 = t327 * t238;
t235 = sin(qJ(2));
t323 = t327 * t235;
t322 = Icges(4,4) - Icges(3,5);
t321 = Icges(4,5) - Icges(3,6);
t320 = t325 * t235 + t324;
t319 = t326 * t238 - t323;
t318 = Icges(4,1) + Icges(3,3);
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t317 = t320 * t236 + t321 * t239;
t316 = -t321 * t236 + t320 * t239;
t315 = t319 * t236 + t322 * t239;
t314 = -t322 * t236 + t319 * t239;
t313 = t325 * t238 - t323;
t312 = t326 * t235 + t324;
t311 = t321 * t235 - t322 * t238;
t310 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t309 = Icges(5,4) - Icges(7,4) - Icges(6,5);
t308 = -Icges(7,5) + Icges(6,4) + Icges(5,5);
t307 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t306 = Icges(6,6) - Icges(7,6) - Icges(5,6);
t305 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t304 = rSges(7,1) + pkin(5);
t303 = -rSges(7,3) - qJ(6);
t225 = -qJD(2) * t239 + V_base(5);
t226 = qJD(2) * t236 + V_base(4);
t230 = V_base(6) + qJD(1);
t302 = (t235 * t313 + t238 * t312) * t230 + (-t235 * t316 + t238 * t314) * t226 + (-t235 * t317 + t238 * t315) * t225;
t301 = (-t322 * t235 - t321 * t238) * t230 + (t236 * t318 + t311 * t239) * t226 + (t311 * t236 - t239 * t318) * t225;
t237 = cos(qJ(4));
t275 = t237 * t239;
t234 = sin(qJ(4));
t279 = t234 * t236;
t188 = -t235 * t275 + t279;
t277 = t236 * t237;
t278 = t234 * t239;
t189 = t235 * t278 + t277;
t274 = t238 * t239;
t300 = t188 * t306 + t189 * t308 + t274 * t305;
t190 = t235 * t277 + t278;
t191 = t235 * t279 - t275;
t276 = t236 * t238;
t299 = -t190 * t306 + t191 * t308 + t276 * t305;
t298 = t188 * t307 - t189 * t309 + t274 * t306;
t297 = -t190 * t307 - t191 * t309 + t276 * t306;
t296 = -t188 * t309 + t189 * t310 + t274 * t308;
t295 = t190 * t309 + t191 * t310 + t276 * t308;
t294 = (-t234 * t308 + t237 * t306) * t238 + t305 * t235;
t293 = (t234 * t309 + t237 * t307) * t238 + t306 * t235;
t292 = (-t234 * t310 - t237 * t309) * t238 + t308 * t235;
t285 = pkin(8) * t235;
t284 = Icges(2,4) * t236;
t273 = t188 * rSges(7,2) + t304 * t189 + t303 * t274;
t272 = -rSges(7,2) * t190 + t304 * t191 + t303 * t276;
t271 = (rSges(7,2) * t237 - t304 * t234) * t238 + t303 * t235;
t260 = pkin(2) * t238 + qJ(3) * t235;
t192 = t260 * t236;
t222 = t236 * pkin(1) - pkin(7) * t239;
t270 = -t192 - t222;
t269 = qJD(3) * t235;
t268 = qJD(4) * t238;
t267 = qJD(6) * t238;
t266 = V_base(5) * pkin(6) + V_base(1);
t217 = pkin(2) * t235 - qJ(3) * t238;
t263 = t225 * t217 + t239 * t269 + t266;
t262 = rSges(3,1) * t238 - rSges(3,2) * t235;
t261 = -rSges(4,2) * t238 + rSges(4,3) * t235;
t223 = pkin(1) * t239 + t236 * pkin(7);
t253 = -V_base(4) * pkin(6) + t230 * t223 + V_base(2);
t252 = V_base(4) * t222 - t223 * V_base(5) + V_base(3);
t194 = t260 * t239;
t249 = t230 * t194 + t236 * t269 + t253;
t248 = -qJD(3) * t238 + t226 * t192 + t252;
t200 = -pkin(3) * t239 + pkin(8) * t276;
t247 = t225 * t285 + (-t200 + t270) * t230 + t263;
t186 = t236 * t268 + t225;
t193 = (-pkin(4) * t234 + qJ(5) * t237) * t238;
t246 = qJD(5) * t188 + t186 * t193 + t247;
t199 = t236 * pkin(3) + pkin(8) * t274;
t245 = t230 * t199 + (-t217 - t285) * t226 + t249;
t244 = t226 * t200 + (-t194 - t199) * t225 + t248;
t143 = pkin(4) * t189 + qJ(5) * t188;
t216 = qJD(4) * t235 + t230;
t243 = -qJD(5) * t190 + t216 * t143 + t245;
t144 = pkin(4) * t191 - qJ(5) * t190;
t187 = t239 * t268 + t226;
t242 = qJD(5) * t238 * t237 + t187 * t144 + t244;
t232 = Icges(2,4) * t239;
t221 = rSges(2,1) * t239 - t236 * rSges(2,2);
t220 = t236 * rSges(2,1) + rSges(2,2) * t239;
t219 = rSges(3,1) * t235 + rSges(3,2) * t238;
t218 = -rSges(4,2) * t235 - rSges(4,3) * t238;
t215 = Icges(2,1) * t239 - t284;
t214 = Icges(2,1) * t236 + t232;
t212 = -Icges(2,2) * t236 + t232;
t211 = Icges(2,2) * t239 + t284;
t203 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t202 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t201 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t181 = -rSges(4,1) * t239 + t236 * t261;
t180 = t236 * rSges(4,1) + t239 * t261;
t179 = t236 * rSges(3,3) + t239 * t262;
t178 = rSges(5,3) * t235 + (-rSges(5,1) * t234 - rSges(5,2) * t237) * t238;
t177 = rSges(6,2) * t235 + (-rSges(6,1) * t234 + rSges(6,3) * t237) * t238;
t175 = -rSges(3,3) * t239 + t236 * t262;
t148 = V_base(5) * rSges(2,3) - t220 * t230 + t266;
t147 = t221 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t146 = t220 * V_base(4) - t221 * V_base(5) + V_base(3);
t142 = rSges(5,1) * t191 + rSges(5,2) * t190 + rSges(5,3) * t276;
t141 = rSges(6,1) * t191 + rSges(6,2) * t276 - rSges(6,3) * t190;
t139 = t189 * rSges(5,1) - t188 * rSges(5,2) + rSges(5,3) * t274;
t138 = t189 * rSges(6,1) + rSges(6,2) * t274 + t188 * rSges(6,3);
t116 = t219 * t225 + (-t175 - t222) * t230 + t266;
t115 = t179 * t230 - t219 * t226 + t253;
t114 = t175 * t226 - t179 * t225 + t252;
t113 = t218 * t225 + (-t181 + t270) * t230 + t263;
t112 = t180 * t230 + (-t217 - t218) * t226 + t249;
t111 = t181 * t226 + (-t180 - t194) * t225 + t248;
t110 = -t142 * t216 + t178 * t186 + t247;
t109 = t139 * t216 - t178 * t187 + t245;
t108 = -t139 * t186 + t142 * t187 + t244;
t107 = t177 * t186 + (-t141 - t144) * t216 + t246;
t106 = t138 * t216 + (-t177 - t193) * t187 + t243;
t105 = t141 * t187 + (-t138 - t143) * t186 + t242;
t104 = -t239 * t267 + t271 * t186 + (-t144 - t272) * t216 + t246;
t103 = -t236 * t267 + t273 * t216 + (-t193 - t271) * t187 + t243;
t102 = -qJD(6) * t235 + t272 * t187 + (-t143 - t273) * t186 + t242;
t1 = m(6) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(3) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(1) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(2) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + (t302 * t236 - t301 * t239) * t225 / 0.2e1 + (t301 * t236 + t302 * t239) * t226 / 0.2e1 + ((-t236 * t211 + t214 * t239 + Icges(1,4)) * V_base(5) + (-t236 * t212 + t215 * t239 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t211 * t239 + t236 * t214 + Icges(1,2)) * V_base(5) + (t212 * t239 + t236 * t215 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t190 * t293 + t191 * t292 + t276 * t294) * t216 + (-t190 * t298 + t191 * t296 + t276 * t300) * t187 + (-t297 * t190 + t295 * t191 + t299 * t276) * t186) * t186 / 0.2e1 + ((t188 * t293 + t189 * t292 + t274 * t294) * t216 + (t298 * t188 + t296 * t189 + t300 * t274) * t187 + (t188 * t297 + t189 * t295 + t274 * t299) * t186) * t187 / 0.2e1 + (((-t234 * t292 + t237 * t293) * t216 + (-t234 * t296 + t237 * t298) * t187 + (-t234 * t295 + t237 * t297) * t186) * t238 + (t186 * t299 + t187 * t300 + t216 * t294) * t235) * t216 / 0.2e1 + ((t235 * t314 + t238 * t316) * t226 + (t235 * t315 + t238 * t317) * t225 + (t235 * t312 - t238 * t313 + Icges(2,3)) * t230) * t230 / 0.2e1 + t230 * V_base(4) * (Icges(2,5) * t239 - Icges(2,6) * t236) + t230 * V_base(5) * (Icges(2,5) * t236 + Icges(2,6) * t239) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
