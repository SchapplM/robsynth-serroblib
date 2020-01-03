% Calculate kinetic energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:06
% EndTime: 2019-12-31 18:47:07
% DurationCPUTime: 0.91s
% Computational Cost: add. (454->82), mult. (932->135), div. (0->0), fcn. (1033->6), ass. (0->50)
t348 = -Icges(6,4) - Icges(5,5);
t347 = Icges(5,6) - Icges(6,6);
t322 = sin(qJ(3));
t323 = sin(qJ(1));
t324 = cos(qJ(3));
t325 = cos(qJ(1));
t269 = -t323 * t322 - t325 * t324;
t346 = t269 ^ 2;
t270 = t325 * t322 - t323 * t324;
t345 = t270 ^ 2;
t344 = rSges(6,1) + pkin(4);
t343 = -rSges(6,3) - qJ(5);
t290 = sin(qJ(4));
t291 = cos(qJ(4));
t342 = t348 * t290 - t347 * t291;
t340 = Icges(6,2) + Icges(5,3);
t339 = t347 * t290 + t348 * t291;
t338 = t270 * t269;
t289 = qJD(1) - qJD(3);
t337 = t342 * t289;
t336 = t343 * t290 - t344 * t291;
t335 = -t340 * t269 + t339 * t270;
t334 = t339 * t269 + t340 * t270;
t317 = -t269 * rSges(6,2) + t336 * t270;
t316 = t270 * rSges(6,2) + t336 * t269;
t315 = qJD(4) * (-t290 * rSges(5,1) - t291 * rSges(5,2));
t312 = -qJD(2) * t325 + qJD(1) * (t325 * pkin(1) + t323 * qJ(2));
t311 = -rSges(5,1) * t291 + rSges(5,2) * t290;
t296 = qJD(1) * t325 * pkin(2) + t312;
t295 = t289 * (-t269 * pkin(3) + t270 * pkin(7)) + t296;
t294 = t343 * qJD(4) * t291 + (t344 * qJD(4) - qJD(5)) * t290;
t281 = t323 * pkin(1) - t325 * qJ(2);
t288 = qJD(2) * t323;
t293 = t288 + (-t323 * pkin(2) - t281) * qJD(1);
t283 = t325 * rSges(2,1) - t323 * rSges(2,2);
t282 = t323 * rSges(2,1) + t325 * rSges(2,2);
t268 = qJD(1) * (t325 * rSges(3,1) + t323 * rSges(3,3)) + t312;
t267 = t288 + (-t323 * rSges(3,1) + t325 * rSges(3,3) - t281) * qJD(1);
t266 = -t270 * pkin(3) - t269 * pkin(7);
t262 = t270 * rSges(5,3) + t311 * t269;
t260 = -t269 * rSges(5,3) + t311 * t270;
t246 = t289 * (-t269 * rSges(4,1) - t270 * rSges(4,2)) + t296;
t245 = -t289 * (-t270 * rSges(4,1) + t269 * rSges(4,2)) + t293;
t244 = (t260 * t270 + t262 * t269) * qJD(4);
t243 = t289 * t262 - t270 * t315 + t295;
t242 = -t269 * t315 + (-t260 - t266) * t289 + t293;
t241 = t294 * t270 + t316 * t289 + t295;
t240 = (-t266 - t317) * t289 + t294 * t269 + t293;
t239 = qJD(5) * t291 + (t316 * t269 + t317 * t270) * qJD(4);
t1 = m(3) * (t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(4) * (t245 ^ 2 + t246 ^ 2) / 0.2e1 + t289 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(6) * (t239 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + (((Icges(5,2) + Icges(6,3)) * t291 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t290 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t291) * t290) * t289 + (t345 + t346) * t342 * qJD(4)) * t289 / 0.2e1 - (-t269 * t337 + (-t334 * t338 + t335 * t346) * qJD(4)) * qJD(4) * t269 / 0.2e1 + (t270 * t337 + (t334 * t345 - t335 * t338) * qJD(4)) * qJD(4) * t270 / 0.2e1 + (m(2) * (t282 ^ 2 + t283 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
