% Calculate kinetic energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:08
% EndTime: 2019-12-31 17:15:09
% DurationCPUTime: 1.14s
% Computational Cost: add. (542->126), mult. (749->202), div. (0->0), fcn. (636->6), ass. (0->84)
t338 = Icges(4,4) + Icges(5,4);
t337 = Icges(4,1) + Icges(5,1);
t336 = Icges(4,2) + Icges(5,2);
t266 = qJ(2) + qJ(3);
t263 = cos(t266);
t335 = t338 * t263;
t262 = sin(t266);
t334 = t338 * t262;
t333 = Icges(4,5) + Icges(5,5);
t332 = Icges(4,6) + Icges(5,6);
t331 = -t336 * t262 + t335;
t330 = t337 * t263 - t334;
t329 = rSges(5,1) + pkin(3);
t328 = Icges(4,3) + Icges(5,3);
t268 = sin(qJ(1));
t270 = cos(qJ(1));
t327 = t331 * t268 - t332 * t270;
t326 = t332 * t268 + t331 * t270;
t325 = t330 * t268 - t333 * t270;
t324 = t333 * t268 + t330 * t270;
t323 = t336 * t263 + t334;
t322 = t337 * t262 + t335;
t321 = -t332 * t262 + t333 * t263;
t320 = rSges(5,3) + qJ(4);
t319 = -rSges(5,2) * t262 + t329 * t263;
t296 = qJD(2) + qJD(3);
t252 = t296 * t268;
t253 = t296 * t270;
t318 = (t327 * t262 - t325 * t263) * t253 + (-t326 * t262 + t324 * t263) * t252 + (-t323 * t262 + t322 * t263) * qJD(1);
t317 = (-t321 * t268 + t328 * t270) * t253 + (t328 * t268 + t321 * t270) * t252 + (t333 * t262 + t332 * t263) * qJD(1);
t269 = cos(qJ(2));
t312 = t269 * pkin(2);
t267 = sin(qJ(2));
t310 = Icges(3,4) * t267;
t309 = Icges(3,4) * t269;
t304 = t319 * t268 - t320 * t270;
t303 = t320 * t268 + t319 * t270;
t216 = -pkin(6) * t270 + t312 * t268;
t217 = pkin(6) * t268 + t312 * t270;
t297 = qJD(2) * t270;
t298 = qJD(2) * t268;
t302 = t216 * t298 + t217 * t297;
t260 = pkin(1) * t268 - pkin(5) * t270;
t301 = -t216 - t260;
t295 = pkin(2) * qJD(2) * t267;
t294 = rSges(5,2) * t263 + t329 * t262;
t293 = t270 * t295;
t292 = rSges(3,1) * t269 - rSges(3,2) * t267;
t291 = rSges(4,1) * t263 - rSges(4,2) * t262;
t289 = Icges(3,1) * t269 - t310;
t286 = -Icges(3,2) * t267 + t309;
t283 = Icges(3,5) * t269 - Icges(3,6) * t267;
t236 = -Icges(3,6) * t270 + t286 * t268;
t238 = -Icges(3,5) * t270 + t289 * t268;
t280 = t236 * t267 - t238 * t269;
t237 = Icges(3,6) * t268 + t286 * t270;
t239 = Icges(3,5) * t268 + t289 * t270;
t279 = -t237 * t267 + t239 * t269;
t255 = Icges(3,2) * t269 + t310;
t256 = Icges(3,1) * t267 + t309;
t278 = -t255 * t267 + t256 * t269;
t251 = qJD(1) * (pkin(1) * t270 + pkin(5) * t268);
t277 = qJD(1) * t217 - t268 * t295 + t251;
t259 = rSges(2,1) * t270 - rSges(2,2) * t268;
t258 = rSges(2,1) * t268 + rSges(2,2) * t270;
t257 = rSges(3,1) * t267 + rSges(3,2) * t269;
t254 = Icges(3,5) * t267 + Icges(3,6) * t269;
t249 = rSges(4,1) * t262 + rSges(4,2) * t263;
t241 = rSges(3,3) * t268 + t292 * t270;
t240 = -rSges(3,3) * t270 + t292 * t268;
t235 = Icges(3,3) * t268 + t283 * t270;
t234 = -Icges(3,3) * t270 + t283 * t268;
t233 = rSges(4,3) * t268 + t291 * t270;
t231 = -rSges(4,3) * t270 + t291 * t268;
t210 = qJD(1) * t241 - t257 * t298 + t251;
t209 = -t257 * t297 + (-t240 - t260) * qJD(1);
t208 = (t240 * t268 + t241 * t270) * qJD(2);
t207 = qJD(1) * t233 - t249 * t252 + t277;
t206 = -t293 - t249 * t253 + (-t231 + t301) * qJD(1);
t205 = t231 * t252 + t233 * t253 + t302;
t204 = t303 * qJD(1) - qJD(4) * t270 - t294 * t252 + t277;
t203 = -t293 + qJD(4) * t268 - t294 * t253 + (t301 - t304) * qJD(1);
t202 = t304 * t252 + t303 * t253 + t302;
t1 = m(3) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + ((t268 * t254 + t278 * t270) * qJD(1) + (t268 ^ 2 * t235 + (t280 * t270 + (-t234 + t279) * t268) * t270) * qJD(2)) * t298 / 0.2e1 - ((-t270 * t254 + t278 * t268) * qJD(1) + (t270 ^ 2 * t234 + (t279 * t268 + (-t235 + t280) * t270) * t268) * qJD(2)) * t297 / 0.2e1 + m(4) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(5) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + (t317 * t268 + t318 * t270) * t252 / 0.2e1 - (t318 * t268 - t317 * t270) * t253 / 0.2e1 + (m(2) * (t258 ^ 2 + t259 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t237 * t269 + t239 * t267) * t268 - (t236 * t269 + t238 * t267) * t270) * qJD(2) - (t325 * t262 + t327 * t263) * t253 + (t324 * t262 + t326 * t263) * t252 + (t269 * t255 + t267 * t256 + t322 * t262 + t323 * t263) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
