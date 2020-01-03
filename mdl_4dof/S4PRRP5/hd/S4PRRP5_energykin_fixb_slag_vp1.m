% Calculate kinetic energy for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:50
% EndTime: 2019-12-31 16:28:51
% DurationCPUTime: 0.99s
% Computational Cost: add. (446->112), mult. (1083->188), div. (0->0), fcn. (1105->6), ass. (0->66)
t321 = Icges(4,1) + Icges(5,1);
t320 = Icges(4,4) + Icges(5,4);
t319 = Icges(4,5) + Icges(5,5);
t318 = Icges(4,2) + Icges(5,2);
t317 = -Icges(5,6) - Icges(4,6);
t316 = -Icges(5,3) - Icges(4,3);
t266 = cos(pkin(6));
t302 = t266 ^ 2;
t265 = sin(pkin(6));
t303 = t265 ^ 2;
t304 = t302 + t303;
t315 = qJD(2) * t304;
t314 = t265 * t266;
t270 = cos(qJ(3));
t268 = sin(qJ(3));
t271 = cos(qJ(2));
t292 = t268 * t271;
t255 = -t265 * t292 - t266 * t270;
t291 = t270 * t271;
t294 = t266 * t268;
t256 = t265 * t291 - t294;
t269 = sin(qJ(2));
t295 = t265 * t269;
t313 = -t317 * t255 + t319 * t256 - t316 * t295;
t257 = t265 * t270 - t266 * t292;
t296 = t265 * t268;
t258 = t266 * t291 + t296;
t293 = t266 * t269;
t312 = -t317 * t257 + t319 * t258 - t316 * t293;
t311 = t318 * t255 + t320 * t256 - t317 * t295;
t310 = t318 * t257 + t320 * t258 - t317 * t293;
t309 = t320 * t255 + t321 * t256 + t319 * t295;
t308 = t320 * t257 + t321 * t258 + t319 * t293;
t307 = t316 * t271 + (t317 * t268 + t319 * t270) * t269;
t306 = t317 * t271 + (-t318 * t268 + t320 * t270) * t269;
t305 = t319 * t271 + (t320 * t268 - t321 * t270) * t269;
t301 = qJD(2) ^ 2;
t298 = pkin(3) * t270;
t273 = qJ(4) * t269 + t271 * t298;
t290 = rSges(5,1) * t256 + rSges(5,2) * t255 + rSges(5,3) * t295 - pkin(3) * t294 + t265 * t273;
t289 = -rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t293 - pkin(3) * t296 - t266 * t273;
t288 = (-qJ(4) - rSges(5,3)) * t271 + (rSges(5,1) * t270 - rSges(5,2) * t268 + t298) * t269;
t287 = qJD(2) * t265;
t286 = qJD(2) * t266;
t285 = qJD(3) * t269;
t284 = qJD(3) * t271;
t283 = qJD(1) + (pkin(2) * t271 + pkin(5) * t269) * t315;
t279 = Icges(3,5) * t271 - Icges(3,6) * t269;
t263 = pkin(2) * t269 - pkin(5) * t271;
t276 = -qJD(2) * t263 + qJD(4) * t269;
t262 = rSges(3,1) * t269 + rSges(3,2) * t271;
t260 = t265 * t285 - t286;
t259 = t266 * t285 + t287;
t254 = -rSges(4,3) * t271 + (rSges(4,1) * t270 - rSges(4,2) * t268) * t269;
t240 = Icges(3,3) * t265 + t266 * t279;
t239 = -Icges(3,3) * t266 + t265 * t279;
t237 = rSges(4,1) * t258 + rSges(4,2) * t257 + rSges(4,3) * t293;
t235 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t295;
t219 = qJD(1) + (rSges(3,1) * t271 - rSges(3,2) * t269) * t315;
t218 = t235 * t284 + t254 * t260 - t263 * t286;
t217 = -t237 * t284 - t254 * t259 - t263 * t287;
t216 = t235 * t259 - t237 * t260 + t283;
t215 = t260 * t288 + t266 * t276 + t284 * t290;
t214 = -t259 * t288 + t265 * t276 + t284 * t289;
t213 = -qJD(4) * t271 + t259 * t290 + t260 * t289 + t283;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t304 * t301 * t262 ^ 2 + t219 ^ 2) / 0.2e1 + t301 * t265 * (-t239 * t314 + t303 * t240) / 0.2e1 - t301 * t266 * (t302 * t239 - t240 * t314) / 0.2e1 + m(4) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + ((-t306 * t257 + t305 * t258 - t307 * t293) * t284 + (t311 * t257 + t309 * t258 + t313 * t293) * t260 + (t310 * t257 + t308 * t258 + t312 * t293) * t259) * t259 / 0.2e1 + ((-t306 * t255 + t305 * t256 - t307 * t295) * t284 + (t311 * t255 + t309 * t256 + t313 * t295) * t260 + (t310 * t255 + t308 * t256 + t312 * t295) * t259) * t260 / 0.2e1 - ((-t312 * t259 - t313 * t260 + t307 * t284) * t271 + ((t306 * t268 + t305 * t270) * t284 + (-t311 * t268 + t309 * t270) * t260 + (-t310 * t268 + t308 * t270) * t259) * t269) * t284 / 0.2e1;
T = t1;
