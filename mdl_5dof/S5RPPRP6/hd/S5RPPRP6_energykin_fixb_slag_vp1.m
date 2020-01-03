% Calculate kinetic energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:57
% EndTime: 2019-12-31 17:54:57
% DurationCPUTime: 0.92s
% Computational Cost: add. (417->108), mult. (546->163), div. (0->0), fcn. (441->6), ass. (0->64)
t382 = -Icges(5,4) + Icges(6,5);
t381 = Icges(5,1) + Icges(6,1);
t380 = Icges(5,2) + Icges(6,3);
t311 = pkin(7) + qJ(4);
t307 = cos(t311);
t379 = t382 * t307;
t306 = sin(t311);
t378 = t382 * t306;
t377 = Icges(6,4) + Icges(5,5);
t376 = Icges(5,6) - Icges(6,6);
t375 = -t380 * t307 + t378;
t374 = t381 * t306 - t379;
t373 = rSges(6,1) + pkin(4);
t372 = rSges(6,3) + qJ(5);
t371 = Icges(6,2) + Icges(5,3);
t315 = sin(qJ(1));
t316 = cos(qJ(1));
t370 = t375 * t315 - t376 * t316;
t369 = -t376 * t315 - t375 * t316;
t368 = t374 * t315 + t377 * t316;
t367 = t377 * t315 - t374 * t316;
t366 = t380 * t306 + t379;
t365 = t381 * t307 + t378;
t364 = t377 * t306 + t376 * t307;
t363 = t373 * t306 - t372 * t307;
t362 = t364 * t315 + t316 * t371;
t361 = t315 * t371 - t364 * t316;
t360 = -t376 * t306 + t377 * t307;
t359 = qJD(1) * t316 * qJ(3) + qJD(3) * t315;
t358 = t306 * t365 - t307 * t366;
t357 = t306 * t367 - t307 * t369;
t356 = -t306 * t368 + t307 * t370;
t355 = (t372 * t306 + t373 * t307) * qJD(4) - qJD(5) * t307;
t312 = sin(pkin(7));
t351 = pkin(3) * t312;
t344 = rSges(6,2) * t316 + t363 * t315;
t343 = rSges(6,2) * t315 - t363 * t316;
t310 = qJD(2) * t315;
t341 = qJD(3) * t316 + t310;
t340 = qJD(4) * t315;
t301 = qJD(1) * (pkin(1) * t316 + qJ(2) * t315);
t336 = -qJD(2) * t316 + t301;
t335 = qJD(1) * (pkin(6) * t316 + t315 * t351) + t301 + t359;
t302 = pkin(1) * t315 - qJ(2) * t316;
t334 = t316 * t351 - t302 + (-pkin(6) - qJ(3)) * t315;
t313 = cos(pkin(7));
t333 = rSges(4,1) * t312 + rSges(4,2) * t313;
t332 = rSges(5,1) * t306 + rSges(5,2) * t307;
t304 = rSges(2,1) * t316 - rSges(2,2) * t315;
t303 = rSges(2,1) * t315 + rSges(2,2) * t316;
t300 = rSges(5,1) * t307 - rSges(5,2) * t306;
t287 = rSges(5,3) * t315 - t332 * t316;
t285 = rSges(5,3) * t316 + t332 * t315;
t271 = qJD(1) * (-rSges(3,2) * t316 + rSges(3,3) * t315) + t336;
t270 = t310 + (rSges(3,2) * t315 + rSges(3,3) * t316 - t302) * qJD(1);
t269 = qJD(1) * (rSges(4,3) * t316 + t333 * t315) + t336 + t359;
t268 = (-t302 + t333 * t316 + (-rSges(4,3) - qJ(3)) * t315) * qJD(1) + t341;
t267 = (-t285 * t315 + t287 * t316) * qJD(4);
t266 = qJD(1) * t285 + (-qJD(4) * t300 - qJD(2)) * t316 + t335;
t265 = t300 * t340 + (-t287 + t334) * qJD(1) + t341;
t264 = qJD(5) * t306 + (-t344 * t315 + t343 * t316) * qJD(4);
t263 = t344 * qJD(1) + (-qJD(2) - t355) * t316 + t335;
t262 = t355 * t315 + (t334 - t343) * qJD(1) + t341;
t1 = m(3) * (t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + (((t306 * t370 + t307 * t368) * t316 + (t306 * t369 + t307 * t367) * t315) * qJD(4) + (t306 * t366 + t307 * t365) * qJD(1)) * qJD(1) / 0.2e1 + ((t361 * t315 ^ 2 + (t356 * t316 + (-t357 + t362) * t315) * t316) * qJD(4) + (t360 * t315 - t358 * t316) * qJD(1)) * t340 / 0.2e1 + ((t362 * t316 ^ 2 + (t357 * t315 + (-t356 + t361) * t316) * t315) * qJD(4) + (t358 * t315 + t360 * t316) * qJD(1)) * qJD(4) * t316 / 0.2e1 + (m(2) * (t303 ^ 2 + t304 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t313 ^ 2 + (-0.2e1 * Icges(4,4) * t313 + Icges(4,2) * t312) * t312) * qJD(1) ^ 2 / 0.2e1;
T = t1;
