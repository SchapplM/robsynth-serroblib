% Calculate kinetic energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:58
% EndTime: 2019-12-05 17:37:58
% DurationCPUTime: 0.57s
% Computational Cost: add. (322->129), mult. (527->210), div. (0->0), fcn. (422->6), ass. (0->74)
t253 = sin(qJ(4));
t287 = pkin(4) * t253;
t254 = sin(qJ(1));
t286 = pkin(6) * t254;
t283 = Icges(5,4) * t253;
t255 = cos(qJ(4));
t282 = Icges(5,4) * t255;
t252 = qJ(4) + qJ(5);
t250 = sin(t252);
t281 = Icges(6,4) * t250;
t251 = cos(t252);
t280 = Icges(6,4) * t251;
t249 = qJD(2) * t254;
t256 = cos(qJ(1));
t279 = qJD(3) * t256 + t249;
t278 = qJD(4) * t254;
t277 = qJD(4) * t256;
t276 = qJD(4) + qJD(5);
t275 = pkin(4) * qJD(4) * t255;
t274 = -qJD(2) * t256 + qJD(1) * (pkin(1) * t256 + qJ(2) * t254);
t273 = rSges(5,1) * t253 + rSges(5,2) * t255;
t272 = rSges(6,1) * t250 + rSges(6,2) * t251;
t271 = Icges(5,1) * t253 + t282;
t270 = Icges(6,1) * t250 + t280;
t269 = Icges(5,2) * t255 + t283;
t268 = Icges(6,2) * t251 + t281;
t267 = Icges(5,5) * t253 + Icges(5,6) * t255;
t266 = Icges(6,5) * t250 + Icges(6,6) * t251;
t224 = Icges(5,6) * t256 + t269 * t254;
t226 = Icges(5,5) * t256 + t271 * t254;
t265 = t224 * t255 + t226 * t253;
t225 = -Icges(5,6) * t254 + t269 * t256;
t227 = -Icges(5,5) * t254 + t271 * t256;
t264 = -t225 * t255 - t227 * t253;
t240 = -Icges(5,2) * t253 + t282;
t241 = Icges(5,1) * t255 - t283;
t263 = t240 * t255 + t241 * t253;
t262 = qJD(1) * t256 * qJ(3) + qJD(3) * t254 + t274;
t242 = pkin(1) * t254 - qJ(2) * t256;
t261 = -pkin(6) * t256 - qJ(3) * t254 - t242;
t237 = t276 * t254;
t238 = t276 * t256;
t260 = qJD(1) * (Icges(6,5) * t251 - Icges(6,6) * t250) + (Icges(6,3) * t256 + t266 * t254) * t238 - (-Icges(6,3) * t254 + t266 * t256) * t237;
t216 = Icges(6,6) * t256 + t268 * t254;
t217 = -Icges(6,6) * t254 + t268 * t256;
t218 = Icges(6,5) * t256 + t270 * t254;
t219 = -Icges(6,5) * t254 + t270 * t256;
t233 = -Icges(6,2) * t250 + t280;
t234 = Icges(6,1) * t251 - t281;
t259 = -(t217 * t251 + t219 * t250) * t237 + (t216 * t251 + t218 * t250) * t238 + (t233 * t251 + t234 * t250) * qJD(1);
t245 = rSges(2,1) * t256 - rSges(2,2) * t254;
t244 = rSges(5,1) * t255 - rSges(5,2) * t253;
t243 = rSges(2,1) * t254 + rSges(2,2) * t256;
t239 = Icges(5,5) * t255 - Icges(5,6) * t253;
t235 = rSges(6,1) * t251 - rSges(6,2) * t250;
t231 = pkin(7) * t256 + t254 * t287;
t230 = -pkin(7) * t254 + t256 * t287;
t229 = -rSges(5,3) * t254 + t273 * t256;
t228 = rSges(5,3) * t256 + t273 * t254;
t223 = -Icges(5,3) * t254 + t267 * t256;
t222 = Icges(5,3) * t256 + t267 * t254;
t221 = -rSges(6,3) * t254 + t272 * t256;
t220 = rSges(6,3) * t256 + t272 * t254;
t213 = qJD(1) * (-rSges(3,2) * t256 + rSges(3,3) * t254) + t274;
t212 = t249 + (rSges(3,2) * t254 + rSges(3,3) * t256 - t242) * qJD(1);
t211 = qJD(1) * (rSges(4,2) * t254 + rSges(4,3) * t256) + t262;
t210 = (t256 * rSges(4,2) - t242 + (-rSges(4,3) - qJ(3)) * t254) * qJD(1) + t279;
t209 = (-t228 * t254 - t229 * t256) * qJD(4);
t208 = t244 * t278 + (t229 - t286) * qJD(1) + t262;
t207 = t244 * t277 + (-t228 + t261) * qJD(1) + t279;
t206 = t254 * t275 + t235 * t237 + (t221 + t230 - t286) * qJD(1) + t262;
t205 = t256 * t275 + t235 * t238 + (-t220 - t231 + t261) * qJD(1) + t279;
t204 = -t220 * t237 - t221 * t238 + (-t230 * t256 - t231 * t254) * qJD(4);
t1 = m(3) * (t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(4) * (t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(5) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 - ((-t254 * t239 + t263 * t256) * qJD(1) + (t254 ^ 2 * t223 + (t265 * t256 + (-t222 + t264) * t254) * t256) * qJD(4)) * t278 / 0.2e1 + ((t256 * t239 + t263 * t254) * qJD(1) + (t256 ^ 2 * t222 + (t264 * t254 + (-t223 + t265) * t256) * t254) * qJD(4)) * t277 / 0.2e1 + m(6) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 - t237 * (-t260 * t254 + t259 * t256) / 0.2e1 + t238 * (t259 * t254 + t260 * t256) / 0.2e1 + ((-(-t253 * t225 + t255 * t227) * t254 + (-t224 * t253 + t226 * t255) * t256) * qJD(4) - (-t217 * t250 + t219 * t251) * t237 + (-t216 * t250 + t218 * t251) * t238 + (-t250 * t233 + t251 * t234 - t253 * t240 + t255 * t241) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t243 ^ 2 + t245 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
