% Calculate kinetic energy for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:52
% EndTime: 2019-12-31 17:12:53
% DurationCPUTime: 0.73s
% Computational Cost: add. (446->85), mult. (474->138), div. (0->0), fcn. (386->6), ass. (0->59)
t288 = Icges(4,4) + Icges(5,4);
t287 = Icges(4,1) + Icges(5,1);
t286 = Icges(4,2) + Icges(5,2);
t226 = cos(qJ(3));
t285 = t288 * t226;
t224 = sin(qJ(3));
t284 = t288 * t224;
t283 = Icges(4,5) + Icges(5,5);
t282 = Icges(4,6) + Icges(5,6);
t281 = -t286 * t224 + t285;
t280 = t287 * t226 - t284;
t279 = rSges(5,1) + pkin(3);
t278 = Icges(4,3) + Icges(5,3);
t222 = qJ(1) + qJ(2);
t219 = sin(t222);
t220 = cos(t222);
t277 = t281 * t219 - t282 * t220;
t276 = t282 * t219 + t281 * t220;
t275 = -t280 * t219 + t283 * t220;
t274 = t283 * t219 + t280 * t220;
t273 = t286 * t226 + t284;
t272 = t287 * t224 + t285;
t271 = -t282 * t224 + t283 * t226;
t270 = rSges(5,3) + qJ(4);
t269 = -rSges(5,2) * t224 + t279 * t226;
t268 = t271 * t219 - t278 * t220;
t267 = t278 * t219 + t271 * t220;
t266 = t283 * t224 + t282 * t226;
t265 = -t273 * t224 + t272 * t226;
t264 = -t276 * t224 + t274 * t226;
t263 = t277 * t224 + t275 * t226;
t256 = pkin(1) * qJD(1);
t251 = t269 * t219 - t270 * t220;
t250 = t270 * t219 + t269 * t220;
t227 = cos(qJ(1));
t218 = t227 * t256;
t221 = qJD(1) + qJD(2);
t249 = t221 * (pkin(2) * t220 + pkin(6) * t219) + t218;
t248 = qJD(3) * t219;
t247 = qJD(3) * t220;
t225 = sin(qJ(1));
t246 = t225 * t256;
t243 = rSges(4,1) * t226 - rSges(4,2) * t224;
t241 = qJD(3) * (-rSges(5,2) * t226 - t279 * t224);
t216 = rSges(2,1) * t227 - rSges(2,2) * t225;
t215 = rSges(2,1) * t225 + rSges(2,2) * t227;
t214 = rSges(4,1) * t224 + rSges(4,2) * t226;
t206 = pkin(2) * t219 - pkin(6) * t220;
t204 = t218 + t221 * (rSges(3,1) * t220 - rSges(3,2) * t219);
t203 = -t246 - t221 * (rSges(3,1) * t219 + rSges(3,2) * t220);
t202 = rSges(4,3) * t219 + t243 * t220;
t200 = -rSges(4,3) * t220 + t243 * t219;
t184 = (t200 * t219 + t202 * t220) * qJD(3);
t183 = t202 * t221 - t214 * t248 + t249;
t182 = -t246 - t214 * t247 + (-t200 - t206) * t221;
t181 = -qJD(4) * t220 + t219 * t241 + t250 * t221 + t249;
t180 = -t246 + qJD(4) * t219 + t220 * t241 + (-t206 - t251) * t221;
t179 = (t251 * t219 + t250 * t220) * qJD(3);
t1 = m(3) * (t203 ^ 2 + t204 ^ 2) / 0.2e1 + t221 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t272 * t224 + t273 * t226) * t221 + ((t275 * t224 - t277 * t226) * t220 + (t274 * t224 + t276 * t226) * t219) * qJD(3)) * t221 / 0.2e1 + (m(2) * (t215 ^ 2 + t216 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t266 * t219 + t265 * t220) * t221 + (t267 * t219 ^ 2 + (t263 * t220 + (t264 - t268) * t219) * t220) * qJD(3)) * t248 / 0.2e1 - ((t265 * t219 - t266 * t220) * t221 + (t268 * t220 ^ 2 + (t264 * t219 + (t263 - t267) * t220) * t219) * qJD(3)) * t247 / 0.2e1;
T = t1;
