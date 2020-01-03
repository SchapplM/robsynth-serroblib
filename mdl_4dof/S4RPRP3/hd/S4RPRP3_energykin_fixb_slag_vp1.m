% Calculate kinetic energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:30
% EndTime: 2019-12-31 16:42:31
% DurationCPUTime: 0.82s
% Computational Cost: add. (428->86), mult. (475->135), div. (0->0), fcn. (386->6), ass. (0->57)
t289 = Icges(4,4) + Icges(5,4);
t288 = Icges(4,1) + Icges(5,1);
t287 = Icges(4,2) + Icges(5,2);
t227 = cos(qJ(3));
t286 = t289 * t227;
t225 = sin(qJ(3));
t285 = t289 * t225;
t284 = Icges(4,5) + Icges(5,5);
t283 = Icges(4,6) + Icges(5,6);
t282 = -t287 * t225 + t286;
t281 = t288 * t227 - t285;
t280 = rSges(5,1) + pkin(3);
t279 = Icges(4,3) + Icges(5,3);
t223 = qJ(1) + pkin(6);
t221 = sin(t223);
t222 = cos(t223);
t278 = t282 * t221 - t283 * t222;
t277 = t283 * t221 + t282 * t222;
t276 = -t281 * t221 + t284 * t222;
t275 = t284 * t221 + t281 * t222;
t274 = t287 * t227 + t285;
t273 = t288 * t225 + t286;
t272 = -t283 * t225 + t284 * t227;
t271 = rSges(5,3) + qJ(4);
t270 = -rSges(5,2) * t225 + t280 * t227;
t269 = t272 * t221 - t279 * t222;
t268 = t279 * t221 + t272 * t222;
t267 = t284 * t225 + t283 * t227;
t266 = -t274 * t225 + t273 * t227;
t265 = -t277 * t225 + t275 * t227;
t264 = t278 * t225 + t276 * t227;
t226 = sin(qJ(1));
t260 = pkin(1) * t226;
t252 = t270 * t221 - t271 * t222;
t251 = t271 * t221 + t270 * t222;
t228 = cos(qJ(1));
t220 = qJD(1) * t228 * pkin(1);
t250 = qJD(1) * (pkin(2) * t222 + pkin(5) * t221) + t220;
t249 = qJD(3) * t221;
t248 = qJD(3) * t222;
t245 = -pkin(2) * t221 + pkin(5) * t222 - t260;
t244 = rSges(4,1) * t227 - rSges(4,2) * t225;
t242 = qJD(3) * (-rSges(5,2) * t227 - t280 * t225);
t218 = rSges(2,1) * t228 - rSges(2,2) * t226;
t217 = rSges(2,1) * t226 + rSges(2,2) * t228;
t216 = rSges(4,1) * t225 + rSges(4,2) * t227;
t206 = t220 + qJD(1) * (rSges(3,1) * t222 - rSges(3,2) * t221);
t205 = (-rSges(3,1) * t221 - rSges(3,2) * t222 - t260) * qJD(1);
t204 = rSges(4,3) * t221 + t244 * t222;
t202 = -rSges(4,3) * t222 + t244 * t221;
t186 = qJD(1) * t204 - t216 * t249 + t250;
t185 = -t216 * t248 + (-t202 + t245) * qJD(1);
t184 = qJD(2) + (t202 * t221 + t204 * t222) * qJD(3);
t183 = t251 * qJD(1) - qJD(4) * t222 + t221 * t242 + t250;
t182 = qJD(4) * t221 + t222 * t242 + (t245 - t252) * qJD(1);
t181 = qJD(2) + (t252 * t221 + t251 * t222) * qJD(3);
t1 = m(3) * (qJD(2) ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(4) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(5) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + (((t276 * t225 - t278 * t227) * t222 + (t275 * t225 + t277 * t227) * t221) * qJD(3) + (t273 * t225 + t274 * t227) * qJD(1)) * qJD(1) / 0.2e1 + ((t268 * t221 ^ 2 + (t264 * t222 + (t265 - t269) * t221) * t222) * qJD(3) + (t221 * t267 + t222 * t266) * qJD(1)) * t249 / 0.2e1 - ((t269 * t222 ^ 2 + (t265 * t221 + (t264 - t268) * t222) * t221) * qJD(3) + (t221 * t266 - t222 * t267) * qJD(1)) * t248 / 0.2e1 + (m(2) * (t217 ^ 2 + t218 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
