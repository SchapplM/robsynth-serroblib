% Calculate kinetic energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:41
% EndTime: 2022-01-20 11:48:42
% DurationCPUTime: 1.27s
% Computational Cost: add. (964->138), mult. (778->218), div. (0->0), fcn. (648->8), ass. (0->93)
t389 = Icges(5,4) + Icges(6,4);
t388 = Icges(5,1) + Icges(6,1);
t387 = Icges(5,2) + Icges(6,2);
t313 = qJ(3) + qJ(4);
t308 = cos(t313);
t386 = t389 * t308;
t306 = sin(t313);
t385 = t389 * t306;
t384 = Icges(5,5) + Icges(6,5);
t383 = Icges(5,6) + Icges(6,6);
t382 = -t387 * t306 + t386;
t381 = t388 * t308 - t385;
t380 = rSges(6,1) + pkin(4);
t379 = Icges(5,3) + Icges(6,3);
t314 = qJ(1) + qJ(2);
t307 = sin(t314);
t309 = cos(t314);
t378 = t382 * t307 - t383 * t309;
t377 = t383 * t307 + t382 * t309;
t376 = t381 * t307 - t384 * t309;
t375 = t384 * t307 + t381 * t309;
t374 = t387 * t308 + t385;
t373 = t388 * t306 + t386;
t372 = -t383 * t306 + t384 * t308;
t371 = rSges(6,3) + qJ(5);
t370 = -rSges(6,2) * t306 + t380 * t308;
t345 = qJD(3) + qJD(4);
t286 = t345 * t307;
t287 = t345 * t309;
t312 = qJD(1) + qJD(2);
t369 = (-t374 * t306 + t373 * t308) * t312 + (t378 * t306 - t376 * t308) * t287 + (-t377 * t306 + t375 * t308) * t286;
t368 = (t384 * t306 + t383 * t308) * t312 + (-t372 * t307 + t379 * t309) * t287 + (t379 * t307 + t372 * t309) * t286;
t317 = cos(qJ(3));
t363 = t317 * pkin(3);
t361 = pkin(1) * qJD(1);
t315 = sin(qJ(3));
t360 = Icges(4,4) * t315;
t359 = Icges(4,4) * t317;
t354 = t370 * t307 - t371 * t309;
t353 = t371 * t307 + t370 * t309;
t257 = -pkin(8) * t309 + t363 * t307;
t258 = pkin(8) * t307 + t363 * t309;
t346 = qJD(3) * t309;
t347 = qJD(3) * t307;
t352 = t257 * t347 + t258 * t346;
t296 = t307 * pkin(2) - t309 * pkin(7);
t351 = -t257 - t296;
t318 = cos(qJ(1));
t305 = t318 * t361;
t350 = t312 * (t309 * pkin(2) + t307 * pkin(7)) + t305;
t344 = pkin(3) * qJD(3) * t315;
t316 = sin(qJ(1));
t343 = t316 * t361;
t342 = t308 * rSges(6,2) + t380 * t306;
t341 = rSges(4,1) * t317 - rSges(4,2) * t315;
t340 = rSges(5,1) * t308 - rSges(5,2) * t306;
t338 = Icges(4,1) * t317 - t360;
t335 = -Icges(4,2) * t315 + t359;
t332 = Icges(4,5) * t317 - Icges(4,6) * t315;
t277 = -Icges(4,6) * t309 + t335 * t307;
t279 = -Icges(4,5) * t309 + t338 * t307;
t329 = t277 * t315 - t279 * t317;
t278 = Icges(4,6) * t307 + t335 * t309;
t280 = Icges(4,5) * t307 + t338 * t309;
t328 = -t278 * t315 + t280 * t317;
t299 = Icges(4,2) * t317 + t360;
t300 = Icges(4,1) * t315 + t359;
t327 = -t299 * t315 + t300 * t317;
t326 = t312 * t258 - t307 * t344 + t350;
t323 = -t309 * t344 - t343;
t303 = t318 * rSges(2,1) - t316 * rSges(2,2);
t302 = t316 * rSges(2,1) + t318 * rSges(2,2);
t301 = t315 * rSges(4,1) + t317 * rSges(4,2);
t298 = Icges(4,5) * t315 + Icges(4,6) * t317;
t295 = t306 * rSges(5,1) + t308 * rSges(5,2);
t284 = t305 + t312 * (t309 * rSges(3,1) - t307 * rSges(3,2));
t283 = -t343 - t312 * (t307 * rSges(3,1) + t309 * rSges(3,2));
t282 = t307 * rSges(4,3) + t341 * t309;
t281 = -t309 * rSges(4,3) + t341 * t307;
t276 = Icges(4,3) * t307 + t332 * t309;
t275 = -Icges(4,3) * t309 + t332 * t307;
t274 = t307 * rSges(5,3) + t340 * t309;
t272 = -t309 * rSges(5,3) + t340 * t307;
t251 = (t281 * t307 + t282 * t309) * qJD(3);
t250 = t312 * t282 - t301 * t347 + t350;
t249 = -t343 - t301 * t346 + (-t281 - t296) * t312;
t248 = t312 * t274 - t286 * t295 + t326;
t247 = -t287 * t295 + (-t272 + t351) * t312 + t323;
t246 = t286 * t272 + t287 * t274 + t352;
t245 = -qJD(5) * t309 - t342 * t286 + t353 * t312 + t326;
t244 = qJD(5) * t307 - t342 * t287 + (t351 - t354) * t312 + t323;
t243 = t354 * t286 + t353 * t287 + t352;
t1 = m(3) * (t283 ^ 2 + t284 ^ 2) / 0.2e1 + t312 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + ((t307 * t298 + t327 * t309) * t312 + (t307 ^ 2 * t276 + (t329 * t309 + (-t275 + t328) * t307) * t309) * qJD(3)) * t347 / 0.2e1 - ((-t309 * t298 + t327 * t307) * t312 + (t309 ^ 2 * t275 + (t328 * t307 + (-t276 + t329) * t309) * t307) * qJD(3)) * t346 / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + (t368 * t307 + t369 * t309) * t286 / 0.2e1 - (t369 * t307 - t368 * t309) * t287 / 0.2e1 + (m(2) * (t302 ^ 2 + t303 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t317 * t278 + t315 * t280) * t307 - (t317 * t277 + t315 * t279) * t309) * qJD(3) - (t376 * t306 + t378 * t308) * t287 + (t375 * t306 + t377 * t308) * t286 + (t317 * t299 + t315 * t300 + t373 * t306 + t374 * t308) * t312) * t312 / 0.2e1;
T = t1;
