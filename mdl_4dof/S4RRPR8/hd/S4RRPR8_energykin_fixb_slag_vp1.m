% Calculate kinetic energy for
% S4RRPR8
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:46
% EndTime: 2019-12-31 17:07:48
% DurationCPUTime: 1.36s
% Computational Cost: add. (386->145), mult. (929->234), div. (0->0), fcn. (894->6), ass. (0->87)
t331 = Icges(3,4) - Icges(4,5);
t330 = Icges(3,1) + Icges(4,1);
t329 = Icges(3,2) + Icges(4,3);
t270 = cos(qJ(2));
t328 = t331 * t270;
t267 = sin(qJ(2));
t327 = t331 * t267;
t326 = Icges(4,4) + Icges(3,5);
t325 = Icges(3,6) - Icges(4,6);
t324 = t329 * t267 - t328;
t323 = t330 * t270 - t327;
t322 = Icges(4,2) + Icges(3,3);
t268 = sin(qJ(1));
t271 = cos(qJ(1));
t321 = t324 * t268 + t325 * t271;
t320 = -t325 * t268 + t324 * t271;
t319 = -t323 * t268 + t326 * t271;
t318 = t326 * t268 + t323 * t271;
t317 = -t329 * t270 - t327;
t316 = t330 * t267 + t328;
t315 = -t325 * t267 + t326 * t270;
t314 = t315 * t268 - t322 * t271;
t313 = t322 * t268 + t315 * t271;
t312 = t326 * t267 + t325 * t270;
t311 = t317 * t267 + t316 * t270;
t310 = t320 * t267 + t318 * t270;
t309 = -t321 * t267 + t319 * t270;
t305 = pkin(3) * t270;
t287 = pkin(2) * t270 + qJ(3) * t267;
t243 = t287 * t268;
t263 = pkin(1) * t268 - pkin(5) * t271;
t299 = -t243 - t263;
t298 = qJD(2) * t268;
t297 = qJD(2) * t271;
t296 = qJD(3) * t267;
t295 = qJD(2) - qJD(4);
t244 = t287 * t271;
t249 = qJD(1) * (pkin(1) * t271 + pkin(5) * t268);
t294 = qJD(1) * t244 + t268 * t296 + t249;
t258 = pkin(2) * t267 - qJ(3) * t270;
t291 = qJD(2) * (-rSges(4,1) * t267 + rSges(4,3) * t270 - t258);
t290 = -qJD(3) * t270 + t243 * t298 + t244 * t297;
t289 = rSges(3,1) * t270 - rSges(3,2) * t267;
t288 = rSges(4,1) * t270 + rSges(4,3) * t267;
t286 = qJD(2) * (-pkin(3) * t267 - t258);
t266 = sin(qJ(4));
t269 = cos(qJ(4));
t246 = -t266 * t270 + t267 * t269;
t273 = t266 * t267 + t269 * t270;
t265 = t271 * t296;
t262 = rSges(2,1) * t271 - rSges(2,2) * t268;
t261 = rSges(2,1) * t268 + rSges(2,2) * t271;
t260 = rSges(3,1) * t267 + rSges(3,2) * t270;
t251 = t295 * t271;
t250 = t295 * t268;
t248 = -pkin(6) * t268 + t271 * t305;
t247 = pkin(6) * t271 + t268 * t305;
t241 = t273 * t271;
t240 = t246 * t271;
t239 = t273 * t268;
t238 = t246 * t268;
t237 = rSges(3,3) * t268 + t271 * t289;
t236 = rSges(4,2) * t268 + t271 * t288;
t235 = -rSges(3,3) * t271 + t268 * t289;
t234 = -rSges(4,2) * t271 + t268 * t288;
t219 = rSges(5,1) * t246 - rSges(5,2) * t273;
t218 = Icges(5,1) * t246 - Icges(5,4) * t273;
t217 = Icges(5,4) * t246 - Icges(5,2) * t273;
t216 = Icges(5,5) * t246 - Icges(5,6) * t273;
t215 = rSges(5,1) * t241 + rSges(5,2) * t240 - rSges(5,3) * t268;
t214 = rSges(5,1) * t239 + rSges(5,2) * t238 + rSges(5,3) * t271;
t213 = Icges(5,1) * t241 + Icges(5,4) * t240 - Icges(5,5) * t268;
t212 = Icges(5,1) * t239 + Icges(5,4) * t238 + Icges(5,5) * t271;
t211 = Icges(5,4) * t241 + Icges(5,2) * t240 - Icges(5,6) * t268;
t210 = Icges(5,4) * t239 + Icges(5,2) * t238 + Icges(5,6) * t271;
t209 = Icges(5,5) * t241 + Icges(5,6) * t240 - Icges(5,3) * t268;
t208 = Icges(5,5) * t239 + Icges(5,6) * t238 + Icges(5,3) * t271;
t207 = qJD(1) * t237 - t260 * t298 + t249;
t206 = -t260 * t297 + (-t235 - t263) * qJD(1);
t205 = (t235 * t268 + t237 * t271) * qJD(2);
t204 = qJD(1) * t236 + t268 * t291 + t294;
t203 = t265 + t271 * t291 + (-t234 + t299) * qJD(1);
t202 = (t234 * t268 + t236 * t271) * qJD(2) + t290;
t201 = -t219 * t250 + t268 * t286 + (t215 + t248) * qJD(1) + t294;
t200 = -t219 * t251 + t265 + t271 * t286 + (-t214 - t247 + t299) * qJD(1);
t199 = t214 * t250 + t215 * t251 + (t247 * t268 + t248 * t271) * qJD(2) + t290;
t1 = m(3) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(5) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + t250 * ((-t209 * t268 + t211 * t240 + t213 * t241) * t250 - (-t208 * t268 + t210 * t240 + t212 * t241) * t251 + (-t216 * t268 + t217 * t240 + t218 * t241) * qJD(1)) / 0.2e1 - t251 * ((t209 * t271 + t211 * t238 + t213 * t239) * t250 - (t208 * t271 + t210 * t238 + t212 * t239) * t251 + (t216 * t271 + t217 * t238 + t218 * t239) * qJD(1)) / 0.2e1 + (m(2) * (t261 ^ 2 + t262 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t313 * t268 ^ 2 + (t309 * t271 + (t310 - t314) * t268) * t271) * qJD(2) + (t312 * t268 + t311 * t271) * qJD(1)) * t298 / 0.2e1 - ((t314 * t271 ^ 2 + (t310 * t268 + (t309 - t313) * t271) * t268) * qJD(2) + (t311 * t268 - t312 * t271) * qJD(1)) * t297 / 0.2e1 + ((-t211 * t273 + t213 * t246) * t250 - (-t210 * t273 + t212 * t246) * t251 + ((t319 * t267 + t321 * t270) * t271 + (t318 * t267 - t320 * t270) * t268) * qJD(2) + (-t273 * t217 + t246 * t218 + t316 * t267 - t317 * t270) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
