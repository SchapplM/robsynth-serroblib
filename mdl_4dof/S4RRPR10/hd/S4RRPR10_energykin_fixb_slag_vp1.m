% Calculate kinetic energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:08
% EndTime: 2019-12-31 17:11:09
% DurationCPUTime: 1.46s
% Computational Cost: add. (389->156), mult. (933->255), div. (0->0), fcn. (876->6), ass. (0->91)
t341 = Icges(3,4) + Icges(4,6);
t340 = Icges(3,1) + Icges(4,2);
t339 = -Icges(3,2) - Icges(4,3);
t276 = cos(qJ(2));
t338 = t341 * t276;
t273 = sin(qJ(2));
t337 = t341 * t273;
t336 = -Icges(4,4) + Icges(3,5);
t335 = Icges(4,5) - Icges(3,6);
t334 = t339 * t273 + t338;
t333 = -t340 * t276 + t337;
t332 = Icges(4,1) + Icges(3,3);
t274 = sin(qJ(1));
t277 = cos(qJ(1));
t331 = t334 * t274 + t335 * t277;
t330 = -t335 * t274 + t334 * t277;
t329 = t333 * t274 + t336 * t277;
t328 = t336 * t274 - t333 * t277;
t327 = t339 * t276 - t337;
t326 = t340 * t273 + t338;
t325 = t335 * t273 + t336 * t276;
t324 = t325 * t274 - t332 * t277;
t323 = t332 * t274 + t325 * t277;
t322 = t336 * t273 - t335 * t276;
t321 = t327 * t273 + t326 * t276;
t320 = -t330 * t273 + t328 * t276;
t319 = t331 * t273 + t329 * t276;
t272 = sin(qJ(4));
t310 = t272 * t274;
t309 = t272 * t277;
t275 = cos(qJ(4));
t308 = t274 * t275;
t307 = t275 * t277;
t306 = t276 * t274;
t305 = t276 * t277;
t292 = pkin(2) * t276 + qJ(3) * t273;
t250 = t292 * t274;
t268 = pkin(1) * t274 - pkin(5) * t277;
t304 = -t250 - t268;
t303 = qJD(2) * t274;
t302 = qJD(2) * t277;
t301 = qJD(3) * t273;
t300 = qJD(4) * t276;
t251 = t292 * t277;
t256 = qJD(1) * (pkin(1) * t277 + pkin(5) * t274);
t299 = qJD(1) * t251 + t274 * t301 + t256;
t263 = pkin(2) * t273 - qJ(3) * t276;
t296 = qJD(2) * (rSges(4,2) * t273 + rSges(4,3) * t276 - t263);
t295 = -qJD(3) * t276 + t250 * t303 + t251 * t302;
t294 = rSges(3,1) * t276 - rSges(3,2) * t273;
t293 = -rSges(4,2) * t276 + rSges(4,3) * t273;
t291 = qJD(2) * (-pkin(6) * t273 - t263);
t271 = qJD(4) * t273 + qJD(1);
t270 = t277 * t301;
t267 = rSges(2,1) * t277 - rSges(2,2) * t274;
t266 = rSges(2,1) * t274 + rSges(2,2) * t277;
t265 = rSges(3,1) * t273 + rSges(3,2) * t276;
t255 = -pkin(3) * t277 + pkin(6) * t306;
t254 = pkin(3) * t274 + pkin(6) * t305;
t253 = t274 * t300 - t302;
t252 = t277 * t300 + t303;
t249 = t273 * t310 - t307;
t248 = t273 * t308 + t309;
t247 = t273 * t309 + t308;
t246 = t273 * t307 - t310;
t244 = -rSges(4,1) * t277 + t293 * t274;
t243 = rSges(4,1) * t274 + t293 * t277;
t242 = rSges(3,3) * t274 + t294 * t277;
t241 = rSges(5,3) * t273 + (-rSges(5,1) * t272 - rSges(5,2) * t275) * t276;
t240 = -rSges(3,3) * t277 + t294 * t274;
t229 = Icges(5,5) * t273 + (-Icges(5,1) * t272 - Icges(5,4) * t275) * t276;
t226 = Icges(5,6) * t273 + (-Icges(5,4) * t272 - Icges(5,2) * t275) * t276;
t223 = Icges(5,3) * t273 + (-Icges(5,5) * t272 - Icges(5,6) * t275) * t276;
t222 = rSges(5,1) * t249 + rSges(5,2) * t248 + rSges(5,3) * t306;
t221 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t305;
t220 = Icges(5,1) * t249 + Icges(5,4) * t248 + Icges(5,5) * t306;
t219 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t305;
t218 = Icges(5,4) * t249 + Icges(5,2) * t248 + Icges(5,6) * t306;
t217 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t305;
t216 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t306;
t215 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t305;
t214 = qJD(1) * t242 - t265 * t303 + t256;
t213 = -t265 * t302 + (-t240 - t268) * qJD(1);
t212 = (t240 * t274 + t242 * t277) * qJD(2);
t211 = qJD(1) * t243 + t274 * t296 + t299;
t210 = t270 + t277 * t296 + (-t244 + t304) * qJD(1);
t209 = (t243 * t277 + t244 * t274) * qJD(2) + t295;
t208 = qJD(1) * t254 + t221 * t271 - t241 * t252 + t274 * t291 + t299;
t207 = -t222 * t271 + t241 * t253 + t270 + t277 * t291 + (-t255 + t304) * qJD(1);
t206 = -t221 * t253 + t222 * t252 + (t254 * t277 + t255 * t274) * qJD(2) + t295;
t1 = m(3) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(4) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(5) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + t252 * ((t215 * t305 + t246 * t217 + t247 * t219) * t252 + (t216 * t305 + t218 * t246 + t220 * t247) * t253 + (t223 * t305 + t226 * t246 + t229 * t247) * t271) / 0.2e1 + t253 * ((t215 * t306 + t217 * t248 + t219 * t249) * t252 + (t216 * t306 + t248 * t218 + t249 * t220) * t253 + (t223 * t306 + t226 * t248 + t229 * t249) * t271) / 0.2e1 + t271 * ((t215 * t252 + t216 * t253 + t223 * t271) * t273 + ((-t217 * t275 - t219 * t272) * t252 + (-t218 * t275 - t220 * t272) * t253 + (-t226 * t275 - t229 * t272) * t271) * t276) / 0.2e1 + (m(2) * (t266 ^ 2 + t267 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t329 * t273 - t331 * t276) * t277 + (t328 * t273 + t330 * t276) * t274) * qJD(2) + (t326 * t273 - t327 * t276) * qJD(1)) * qJD(1) / 0.2e1 + ((t323 * t274 ^ 2 + (t319 * t277 + (t320 - t324) * t274) * t277) * qJD(2) + (t322 * t274 + t321 * t277) * qJD(1)) * t303 / 0.2e1 - ((t324 * t277 ^ 2 + (t320 * t274 + (t319 - t323) * t277) * t274) * qJD(2) + (t321 * t274 - t322 * t277) * qJD(1)) * t302 / 0.2e1;
T = t1;
