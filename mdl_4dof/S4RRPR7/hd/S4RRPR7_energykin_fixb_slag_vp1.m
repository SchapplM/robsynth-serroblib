% Calculate kinetic energy for
% S4RRPR7
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:56
% EndTime: 2019-12-31 17:05:57
% DurationCPUTime: 1.13s
% Computational Cost: add. (651->175), mult. (937->284), div. (0->0), fcn. (880->8), ass. (0->99)
t323 = Icges(3,3) + Icges(4,3);
t265 = qJ(2) + pkin(7);
t262 = sin(t265);
t263 = cos(t265);
t268 = sin(qJ(2));
t271 = cos(qJ(2));
t322 = Icges(3,5) * t271 + Icges(4,5) * t263 - Icges(3,6) * t268 - Icges(4,6) * t262;
t269 = sin(qJ(1));
t272 = cos(qJ(1));
t321 = t322 * t269 - t323 * t272;
t320 = t323 * t269 + t322 * t272;
t319 = Icges(3,5) * t268 + Icges(4,5) * t262 + Icges(3,6) * t271 + Icges(4,6) * t263;
t306 = Icges(4,4) * t262;
t248 = Icges(4,2) * t263 + t306;
t305 = Icges(4,4) * t263;
t249 = Icges(4,1) * t262 + t305;
t308 = Icges(3,4) * t268;
t254 = Icges(3,2) * t271 + t308;
t307 = Icges(3,4) * t271;
t255 = Icges(3,1) * t268 + t307;
t318 = -t248 * t262 + t249 * t263 - t254 * t268 + t255 * t271;
t282 = -Icges(4,2) * t262 + t305;
t226 = Icges(4,6) * t269 + t282 * t272;
t284 = Icges(4,1) * t263 - t306;
t228 = Icges(4,5) * t269 + t284 * t272;
t283 = -Icges(3,2) * t268 + t307;
t234 = Icges(3,6) * t269 + t283 * t272;
t285 = Icges(3,1) * t271 - t308;
t236 = Icges(3,5) * t269 + t285 * t272;
t317 = -t226 * t262 + t228 * t263 - t234 * t268 + t236 * t271;
t225 = -Icges(4,6) * t272 + t282 * t269;
t227 = -Icges(4,5) * t272 + t284 * t269;
t233 = -Icges(3,6) * t272 + t283 * t269;
t235 = -Icges(3,5) * t272 + t285 * t269;
t316 = t225 * t262 - t227 * t263 + t233 * t268 - t235 * t271;
t312 = pkin(2) * t268;
t310 = t271 * pkin(2);
t304 = t262 * t269;
t303 = t262 * t272;
t267 = sin(qJ(4));
t302 = t269 * t267;
t270 = cos(qJ(4));
t301 = t269 * t270;
t300 = t272 * t267;
t299 = t272 * t270;
t221 = -qJ(3) * t272 + t310 * t269;
t222 = qJ(3) * t269 + t310 * t272;
t295 = qJD(2) * t272;
t296 = qJD(2) * t269;
t298 = t221 * t296 + t222 * t295;
t259 = t269 * pkin(1) - t272 * pkin(5);
t297 = -t221 - t259;
t294 = qJD(4) * t262;
t291 = pkin(3) * t263 + pkin(6) * t262;
t252 = qJD(1) * (t272 * pkin(1) + t269 * pkin(5));
t290 = qJD(1) * t222 - qJD(3) * t272 + t252;
t289 = rSges(3,1) * t271 - rSges(3,2) * t268;
t288 = rSges(4,1) * t263 - rSges(4,2) * t262;
t287 = qJD(2) * (-t262 * rSges(4,1) - t263 * rSges(4,2) - t312);
t286 = qJD(2) * (-t262 * pkin(3) + t263 * pkin(6) - t312);
t264 = qJD(3) * t269;
t260 = -qJD(4) * t263 + qJD(1);
t258 = t272 * rSges(2,1) - t269 * rSges(2,2);
t257 = t269 * rSges(2,1) + t272 * rSges(2,2);
t256 = t268 * rSges(3,1) + t271 * rSges(3,2);
t246 = t269 * t294 - t295;
t245 = t272 * t294 + t296;
t244 = t263 * t299 + t302;
t243 = -t263 * t300 + t301;
t242 = t263 * t301 - t300;
t241 = -t263 * t302 - t299;
t240 = t291 * t272;
t239 = t291 * t269;
t238 = t269 * rSges(3,3) + t289 * t272;
t237 = -t272 * rSges(3,3) + t289 * t269;
t230 = t269 * rSges(4,3) + t288 * t272;
t229 = -t272 * rSges(4,3) + t288 * t269;
t220 = -t263 * rSges(5,3) + (rSges(5,1) * t270 - rSges(5,2) * t267) * t262;
t218 = -Icges(5,5) * t263 + (Icges(5,1) * t270 - Icges(5,4) * t267) * t262;
t217 = -Icges(5,6) * t263 + (Icges(5,4) * t270 - Icges(5,2) * t267) * t262;
t216 = -Icges(5,3) * t263 + (Icges(5,5) * t270 - Icges(5,6) * t267) * t262;
t213 = qJD(1) * t238 - t256 * t296 + t252;
t212 = -t256 * t295 + (-t237 - t259) * qJD(1);
t211 = (t237 * t269 + t238 * t272) * qJD(2);
t210 = t244 * rSges(5,1) + t243 * rSges(5,2) + rSges(5,3) * t303;
t209 = t242 * rSges(5,1) + t241 * rSges(5,2) + rSges(5,3) * t304;
t208 = Icges(5,1) * t244 + Icges(5,4) * t243 + Icges(5,5) * t303;
t207 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t304;
t206 = Icges(5,4) * t244 + Icges(5,2) * t243 + Icges(5,6) * t303;
t205 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t304;
t204 = Icges(5,5) * t244 + Icges(5,6) * t243 + Icges(5,3) * t303;
t203 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t304;
t202 = qJD(1) * t230 + t269 * t287 + t290;
t201 = t264 + t272 * t287 + (-t229 + t297) * qJD(1);
t200 = (t229 * t269 + t230 * t272) * qJD(2) + t298;
t199 = qJD(1) * t240 + t260 * t210 - t245 * t220 + t269 * t286 + t290;
t198 = -t260 * t209 + t246 * t220 + t264 + t272 * t286 + (-t239 + t297) * qJD(1);
t197 = t245 * t209 - t246 * t210 + (t239 * t269 + t240 * t272) * qJD(2) + t298;
t1 = m(3) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(4) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(5) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + t245 * ((t204 * t303 + t243 * t206 + t244 * t208) * t245 + (t203 * t303 + t243 * t205 + t244 * t207) * t246 + (t216 * t303 + t243 * t217 + t244 * t218) * t260) / 0.2e1 + t246 * ((t204 * t304 + t241 * t206 + t242 * t208) * t245 + (t203 * t304 + t241 * t205 + t242 * t207) * t246 + (t216 * t304 + t241 * t217 + t242 * t218) * t260) / 0.2e1 + t260 * ((-t203 * t246 - t204 * t245 - t216 * t260) * t263 + ((-t206 * t267 + t208 * t270) * t245 + (-t205 * t267 + t207 * t270) * t246 + (-t217 * t267 + t218 * t270) * t260) * t262) / 0.2e1 + (m(2) * (t257 ^ 2 + t258 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t263 * t225 - t262 * t227 - t271 * t233 - t268 * t235) * t272 + (t263 * t226 + t262 * t228 + t271 * t234 + t268 * t236) * t269) * qJD(2) + (t263 * t248 + t262 * t249 + t271 * t254 + t268 * t255) * qJD(1)) * qJD(1) / 0.2e1 + ((t320 * t269 ^ 2 + (t316 * t272 + (t317 - t321) * t269) * t272) * qJD(2) + (t319 * t269 + t318 * t272) * qJD(1)) * t296 / 0.2e1 - ((t321 * t272 ^ 2 + (t317 * t269 + (t316 - t320) * t272) * t269) * qJD(2) + (t318 * t269 - t319 * t272) * qJD(1)) * t295 / 0.2e1;
T = t1;
