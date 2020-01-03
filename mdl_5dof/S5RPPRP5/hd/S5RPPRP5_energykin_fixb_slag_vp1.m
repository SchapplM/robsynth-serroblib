% Calculate kinetic energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:19
% EndTime: 2019-12-31 17:53:20
% DurationCPUTime: 0.71s
% Computational Cost: add. (423->124), mult. (1039->186), div. (0->0), fcn. (1103->6), ass. (0->67)
t382 = Icges(5,1) + Icges(6,1);
t381 = Icges(5,4) - Icges(6,5);
t380 = Icges(6,4) + Icges(5,5);
t379 = Icges(5,2) + Icges(6,3);
t378 = Icges(6,2) + Icges(5,3);
t377 = -Icges(5,6) + Icges(6,6);
t376 = rSges(6,1) + pkin(4);
t375 = rSges(6,3) + qJ(5);
t336 = sin(qJ(4));
t337 = sin(qJ(1));
t334 = sin(pkin(7));
t362 = cos(qJ(4));
t347 = t334 * t362;
t335 = cos(pkin(7));
t360 = t335 * t337;
t315 = t336 * t360 - t337 * t347;
t322 = t334 * t336 + t335 * t362;
t316 = t322 * t337;
t338 = cos(qJ(1));
t359 = t335 * t338;
t317 = t336 * t359 - t338 * t347;
t318 = t322 * t338;
t374 = -(t377 * t317 + t380 * t318 - t378 * t337) * t337 + (t377 * t315 + t380 * t316 + t378 * t338) * t338;
t373 = t379 * t315 - t381 * t316 + t377 * t338;
t372 = -t379 * t317 + t381 * t318 + t377 * t337;
t369 = -t381 * t315 + t382 * t316 + t380 * t338;
t368 = t381 * t317 - t382 * t318 + t380 * t337;
t323 = -t335 * t336 + t347;
t367 = t379 * t322 - t381 * t323;
t366 = t377 * t322 + t380 * t323;
t365 = -t381 * t322 + t382 * t323;
t364 = t335 ^ 2;
t358 = -rSges(6,2) * t338 - t375 * t315 - t376 * t316;
t357 = -rSges(6,2) * t337 + t375 * t317 + t376 * t318;
t326 = t337 * pkin(1) - t338 * qJ(2);
t341 = pkin(2) * t335 + qJ(3) * t334;
t356 = -t341 * t337 - t326;
t333 = qJD(2) * t337;
t353 = qJD(3) * t334;
t355 = t338 * t353 + t333;
t354 = qJD(1) * t337;
t352 = qJD(3) * t335;
t351 = qJD(4) * t337;
t350 = qJD(4) * t338;
t325 = qJD(1) * (t338 * pkin(1) + t337 * qJ(2));
t349 = qJD(1) * t341 * t338 + t337 * t353 + t325;
t348 = -pkin(3) * t360 - t338 * pkin(6) + t356;
t344 = qJD(4) * (t375 * t322 + t376 * t323);
t343 = rSges(3,1) * t335 - rSges(3,2) * t334;
t342 = rSges(4,1) * t335 + rSges(4,3) * t334;
t340 = -qJD(2) * t338 + qJD(1) * (pkin(3) * t359 - t337 * pkin(6)) + t349;
t328 = t338 * rSges(2,1) - t337 * rSges(2,2);
t327 = t337 * rSges(2,1) + t338 * rSges(2,2);
t314 = t323 * rSges(5,1) - t322 * rSges(5,2);
t305 = rSges(3,3) * t354 + t325 + (qJD(1) * t343 - qJD(2)) * t338;
t304 = t333 + (t338 * rSges(3,3) - t343 * t337 - t326) * qJD(1);
t301 = rSges(5,1) * t318 - rSges(5,2) * t317 - rSges(5,3) * t337;
t299 = rSges(5,1) * t316 - rSges(5,2) * t315 + rSges(5,3) * t338;
t285 = rSges(4,2) * t354 + (qJD(1) * t342 - qJD(2)) * t338 + t349;
t284 = (t338 * rSges(4,2) - t342 * t337 + t356) * qJD(1) + t355;
t283 = -t352 + (-t299 * t337 - t301 * t338) * qJD(4);
t282 = qJD(1) * t301 + t314 * t351 + t340;
t281 = t314 * t350 + (-t299 + t348) * qJD(1) + t355;
t280 = -t352 + qJD(5) * t322 + (t358 * t337 - t357 * t338) * qJD(4);
t279 = t357 * qJD(1) + qJD(5) * t315 + t337 * t344 + t340;
t278 = qJD(5) * t317 + t338 * t344 + (t348 + t358) * qJD(1) + t355;
t1 = m(3) * (t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t364 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(5) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(6) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + (((t373 * t322 + t369 * t323) * t338 + (t372 * t322 + t368 * t323) * t337) * qJD(4) + (t367 * t322 + t365 * t323) * qJD(1)) * qJD(1) / 0.2e1 - (((t373 * t317 + t369 * t318) * t338 + (t372 * t317 + t368 * t318 - t374) * t337) * qJD(4) + (t367 * t317 + t365 * t318 - t366 * t337) * qJD(1)) * t351 / 0.2e1 + (((t372 * t315 + t368 * t316) * t337 + (t373 * t315 + t369 * t316 + t374) * t338) * qJD(4) + (t367 * t315 + t365 * t316 + t366 * t338) * qJD(1)) * t350 / 0.2e1 + (m(2) * (t327 ^ 2 + t328 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t364 + ((Icges(3,1) + Icges(4,1)) * t334 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t335) * t334) * qJD(1) ^ 2 / 0.2e1;
T = t1;
