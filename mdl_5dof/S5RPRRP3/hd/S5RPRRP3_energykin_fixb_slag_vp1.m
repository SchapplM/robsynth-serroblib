% Calculate kinetic energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:52
% EndTime: 2020-01-03 11:46:54
% DurationCPUTime: 1.77s
% Computational Cost: add. (935->140), mult. (779->216), div. (0->0), fcn. (648->8), ass. (0->92)
t404 = Icges(5,4) + Icges(6,4);
t403 = Icges(5,1) + Icges(6,1);
t402 = Icges(5,2) + Icges(6,2);
t323 = qJ(3) + qJ(4);
t319 = cos(t323);
t401 = t404 * t319;
t318 = sin(t323);
t400 = t404 * t318;
t399 = Icges(5,5) + Icges(6,5);
t398 = -Icges(5,6) - Icges(6,6);
t397 = -t402 * t318 + t401;
t396 = t403 * t319 - t400;
t395 = rSges(6,1) + pkin(4);
t394 = Icges(5,3) + Icges(6,3);
t322 = qJ(1) + pkin(8);
t316 = sin(t322);
t317 = cos(t322);
t393 = t397 * t316 + t398 * t317;
t392 = t398 * t316 - t397 * t317;
t391 = t396 * t316 - t399 * t317;
t390 = -t399 * t316 - t396 * t317;
t389 = t402 * t319 + t400;
t388 = t403 * t318 + t401;
t387 = t398 * t318 + t399 * t319;
t386 = -rSges(6,3) - qJ(5);
t385 = -rSges(6,2) * t318 + t395 * t319;
t358 = -qJD(3) - qJD(4);
t294 = t358 * t316;
t295 = t358 * t317;
t384 = (t389 * t318 - t388 * t319) * qJD(1) + (t392 * t318 - t390 * t319) * t294 + (t393 * t318 - t391 * t319) * t295;
t383 = (-t387 * t316 + t394 * t317) * t295 + (t394 * t316 + t387 * t317) * t294 + (-t399 * t318 + t398 * t319) * qJD(1);
t326 = cos(qJ(3));
t375 = t326 * pkin(3);
t373 = pkin(1) * qJD(1);
t324 = sin(qJ(3));
t372 = Icges(4,4) * t324;
t371 = Icges(4,4) * t326;
t366 = t385 * t316 + t386 * t317;
t365 = t386 * t316 - t385 * t317;
t266 = -pkin(7) * t316 - t375 * t317;
t296 = -t317 * pkin(2) - t316 * pkin(6);
t364 = -t266 - t296;
t325 = sin(qJ(1));
t314 = t325 * t373;
t363 = qJD(1) * (t316 * pkin(2) - t317 * pkin(6)) + t314;
t360 = qJD(3) * t316;
t359 = qJD(3) * t317;
t357 = pkin(3) * qJD(3) * t324;
t356 = t319 * rSges(6,2) + t395 * t318;
t265 = -pkin(7) * t317 + t375 * t316;
t355 = qJD(1) * t265 + t317 * t357 + t363;
t354 = rSges(4,1) * t326 - rSges(4,2) * t324;
t353 = rSges(5,1) * t319 - rSges(5,2) * t318;
t351 = Icges(4,1) * t326 - t372;
t348 = -Icges(4,2) * t324 + t371;
t345 = Icges(4,5) * t326 - Icges(4,6) * t324;
t285 = -Icges(4,6) * t317 + t348 * t316;
t287 = -Icges(4,5) * t317 + t351 * t316;
t338 = -t285 * t324 + t287 * t326;
t286 = -Icges(4,6) * t316 - t348 * t317;
t288 = -Icges(4,5) * t316 - t351 * t317;
t337 = t286 * t324 - t288 * t326;
t308 = Icges(4,2) * t326 + t372;
t309 = Icges(4,1) * t324 + t371;
t334 = t308 * t324 - t309 * t326;
t327 = cos(qJ(1));
t315 = t327 * t373;
t333 = -t316 * t357 + t315;
t332 = t265 * t360 - t266 * t359 + qJD(2);
t312 = -t327 * rSges(2,1) + t325 * rSges(2,2);
t311 = t325 * rSges(2,1) + t327 * rSges(2,2);
t310 = t324 * rSges(4,1) + t326 * rSges(4,2);
t307 = Icges(4,5) * t324 + Icges(4,6) * t326;
t304 = t318 * rSges(5,1) + t319 * rSges(5,2);
t292 = t315 - qJD(1) * (-t317 * rSges(3,1) + t316 * rSges(3,2));
t291 = t314 + qJD(1) * (t316 * rSges(3,1) + t317 * rSges(3,2));
t290 = -t316 * rSges(4,3) - t354 * t317;
t289 = -t317 * rSges(4,3) + t354 * t316;
t284 = -Icges(4,3) * t316 - t345 * t317;
t283 = -Icges(4,3) * t317 + t345 * t316;
t282 = -t316 * rSges(5,3) - t353 * t317;
t280 = -t317 * rSges(5,3) + t353 * t316;
t260 = -t310 * t360 + t315 + (-t290 - t296) * qJD(1);
t259 = qJD(1) * t289 + t310 * t359 + t363;
t258 = qJD(2) + (t289 * t316 - t290 * t317) * qJD(3);
t257 = t294 * t304 + (-t282 + t364) * qJD(1) + t333;
t256 = qJD(1) * t280 - t295 * t304 + t355;
t255 = -t294 * t280 + t295 * t282 + t332;
t254 = -qJD(5) * t317 + t356 * t294 + (t364 - t365) * qJD(1) + t333;
t253 = t366 * qJD(1) - qJD(5) * t316 - t356 * t295 + t355;
t252 = -t366 * t294 + t365 * t295 + t332;
t1 = m(3) * (qJD(2) ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + m(4) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 - ((-t317 * t307 - t334 * t316) * qJD(1) + (t317 ^ 2 * t283 + (t337 * t316 + (t284 - t338) * t317) * t316) * qJD(3)) * t359 / 0.2e1 - ((-t316 * t307 + t334 * t317) * qJD(1) + (t316 ^ 2 * t284 + (t338 * t317 + (t283 - t337) * t316) * t317) * qJD(3)) * t360 / 0.2e1 + m(5) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + m(6) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + (t383 * t316 + t384 * t317) * t294 / 0.2e1 + (-t384 * t316 + t383 * t317) * t295 / 0.2e1 + (m(2) * (t311 ^ 2 + t312 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + ((-(t326 * t285 + t324 * t287) * t317 - (t326 * t286 + t324 * t288) * t316) * qJD(3) + (t391 * t318 + t393 * t319) * t295 + (t390 * t318 + t392 * t319) * t294 + (t326 * t308 + t324 * t309 + t388 * t318 + t389 * t319) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
