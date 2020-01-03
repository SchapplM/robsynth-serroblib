% Calculate kinetic energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:25
% EndTime: 2019-12-31 17:34:26
% DurationCPUTime: 0.80s
% Computational Cost: add. (427->69), mult. (894->115), div. (0->0), fcn. (1026->6), ass. (0->43)
t339 = -Icges(5,5) - Icges(6,5);
t338 = Icges(5,6) + Icges(6,6);
t281 = sin(pkin(7));
t282 = cos(pkin(7));
t317 = sin(qJ(3));
t318 = cos(qJ(3));
t266 = -t281 * t317 - t282 * t318;
t337 = t266 ^ 2;
t267 = -t281 * t318 + t282 * t317;
t336 = t267 ^ 2;
t335 = rSges(6,1) + pkin(4);
t284 = sin(qJ(4));
t285 = cos(qJ(4));
t334 = t339 * t284 - t338 * t285;
t332 = Icges(5,3) + Icges(6,3);
t331 = t338 * t284 + t339 * t285;
t330 = rSges(6,3) + qJ(5);
t329 = t267 * t266;
t328 = rSges(6,2) * t284 - t335 * t285;
t327 = t334 * qJD(3);
t326 = -t332 * t266 + t331 * t267;
t325 = t331 * t266 + t332 * t267;
t309 = -t330 * t266 + t328 * t267;
t308 = t328 * t266 + t330 * t267;
t280 = qJD(2) * t281;
t307 = qJD(3) * (-t267 * pkin(3) - t266 * pkin(6)) + t280;
t306 = qJD(2) * t282;
t305 = qJD(4) * (-t284 * rSges(5,1) - t285 * rSges(5,2));
t302 = -rSges(5,1) * t285 + rSges(5,2) * t284;
t300 = qJD(4) * (t285 * rSges(6,2) + t335 * t284);
t287 = qJD(1) ^ 2;
t265 = -t266 * pkin(3) + t267 * pkin(6);
t263 = -t306 - qJD(3) * (-t266 * rSges(4,1) - t267 * rSges(4,2));
t262 = t280 + qJD(3) * (-t267 * rSges(4,1) + t266 * rSges(4,2));
t261 = t267 * rSges(5,3) + t302 * t266;
t259 = -t266 * rSges(5,3) + t302 * t267;
t243 = -t267 * t305 - t306 + (-t261 - t265) * qJD(3);
t242 = qJD(3) * t259 - t266 * t305 + t307;
t241 = qJD(1) + (t259 * t267 + t261 * t266) * qJD(4);
t240 = -t306 - qJD(5) * t266 + t267 * t300 + (-t265 - t308) * qJD(3);
t239 = t309 * qJD(3) + qJD(5) * t267 + t266 * t300 + t307;
t238 = qJD(1) + (t308 * t266 + t309 * t267) * qJD(4);
t1 = m(2) * t287 / 0.2e1 + m(3) * (t287 + (t281 ^ 2 + t282 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t262 ^ 2 + t263 ^ 2 + t287) / 0.2e1 + qJD(3) ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(6) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 - (((-Icges(5,2) - Icges(6,2)) * t285 ^ 2 + ((-Icges(5,1) - Icges(6,1)) * t284 + (-2 * Icges(5,4) - 2 * Icges(6,4)) * t285) * t284) * qJD(3) + (t336 + t337) * t334 * qJD(4)) * qJD(3) / 0.2e1 - ((-t325 * t329 + t326 * t337) * qJD(4) + t266 * t327) * qJD(4) * t266 / 0.2e1 + ((t325 * t336 - t326 * t329) * qJD(4) - t267 * t327) * qJD(4) * t267 / 0.2e1;
T = t1;
