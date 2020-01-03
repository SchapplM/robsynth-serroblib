% Calculate kinetic energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:36
% EndTime: 2019-12-31 19:52:36
% DurationCPUTime: 0.83s
% Computational Cost: add. (473->100), mult. (512->154), div. (0->0), fcn. (407->6), ass. (0->63)
t357 = -Icges(5,4) + Icges(6,5);
t356 = Icges(5,1) + Icges(6,1);
t355 = Icges(5,2) + Icges(6,3);
t291 = cos(qJ(4));
t354 = t357 * t291;
t289 = sin(qJ(4));
t353 = t357 * t289;
t352 = Icges(6,4) + Icges(5,5);
t351 = Icges(5,6) - Icges(6,6);
t350 = -t355 * t291 + t353;
t349 = t356 * t289 - t354;
t348 = rSges(6,1) + pkin(4);
t347 = rSges(6,3) + qJ(5);
t346 = Icges(6,2) + Icges(5,3);
t288 = qJ(1) + qJ(2);
t284 = sin(t288);
t285 = cos(t288);
t345 = t350 * t284 - t351 * t285;
t344 = -t351 * t284 - t350 * t285;
t343 = t349 * t284 + t352 * t285;
t342 = t352 * t284 - t349 * t285;
t341 = t355 * t289 + t354;
t340 = t356 * t291 + t353;
t339 = t352 * t289 + t351 * t291;
t338 = t348 * t289 - t347 * t291;
t337 = t339 * t284 + t346 * t285;
t336 = t346 * t284 - t339 * t285;
t335 = -t351 * t289 + t352 * t291;
t334 = t340 * t289 - t341 * t291;
t333 = t342 * t289 - t344 * t291;
t332 = -t343 * t289 + t345 * t291;
t331 = -(t347 * t289 + t348 * t291) * qJD(4) + qJD(5) * t291;
t287 = qJD(1) + qJD(2);
t325 = pkin(7) * t287;
t324 = pkin(1) * qJD(1);
t319 = t285 * rSges(6,2) + t338 * t284;
t318 = t284 * rSges(6,2) - t338 * t285;
t292 = cos(qJ(1));
t283 = t292 * t324;
t317 = t287 * (t285 * pkin(2) + t284 * qJ(3)) + t283;
t315 = qJD(4) * t284;
t290 = sin(qJ(1));
t313 = t290 * t324;
t312 = t285 * t325 + t317;
t309 = qJD(3) * t284 - t313;
t308 = rSges(5,1) * t289 + rSges(5,2) * t291;
t281 = t292 * rSges(2,1) - t290 * rSges(2,2);
t280 = t291 * rSges(5,1) - t289 * rSges(5,2);
t277 = t290 * rSges(2,1) + t292 * rSges(2,2);
t269 = t284 * pkin(2) - t285 * qJ(3);
t265 = t283 + t287 * (t285 * rSges(3,1) - t284 * rSges(3,2));
t264 = -t313 - t287 * (t284 * rSges(3,1) + t285 * rSges(3,2));
t263 = t284 * rSges(5,3) - t308 * t285;
t261 = t285 * rSges(5,3) + t308 * t284;
t247 = -qJD(3) * t285 + t287 * (-t285 * rSges(4,2) + t284 * rSges(4,3)) + t317;
t246 = (t284 * rSges(4,2) + t285 * rSges(4,3) - t269) * t287 + t309;
t245 = (-t261 * t284 + t263 * t285) * qJD(4);
t244 = t287 * t261 + (-qJD(4) * t280 - qJD(3)) * t285 + t312;
t243 = t280 * t315 + (-pkin(7) * t284 - t263 - t269) * t287 + t309;
t242 = qJD(5) * t289 + (-t319 * t284 + t318 * t285) * qJD(4);
t241 = t319 * t287 + (-qJD(3) + t331) * t285 + t312;
t240 = (-t269 - t318) * t287 + (-t325 - t331) * t284 + t309;
t1 = m(3) * (t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(4) * (t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(5) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + m(6) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,1)) * t287 ^ 2 / 0.2e1 + ((t341 * t289 + t340 * t291) * t287 + ((t345 * t289 + t343 * t291) * t285 + (t344 * t289 + t342 * t291) * t284) * qJD(4)) * t287 / 0.2e1 + (m(2) * (t277 ^ 2 + t281 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t335 * t284 - t334 * t285) * t287 + (t336 * t284 ^ 2 + (t332 * t285 + (-t333 + t337) * t284) * t285) * qJD(4)) * t315 / 0.2e1 + ((t334 * t284 + t335 * t285) * t287 + (t337 * t285 ^ 2 + (t333 * t284 + (-t332 + t336) * t285) * t284) * qJD(4)) * qJD(4) * t285 / 0.2e1;
T = t1;
