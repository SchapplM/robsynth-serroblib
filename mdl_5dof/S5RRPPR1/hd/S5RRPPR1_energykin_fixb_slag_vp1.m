% Calculate kinetic energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:09
% EndTime: 2019-12-05 18:18:09
% DurationCPUTime: 0.35s
% Computational Cost: add. (500->96), mult. (318->153), div. (0->0), fcn. (232->10), ass. (0->60)
t248 = qJ(1) + qJ(2);
t243 = sin(t248);
t277 = pkin(2) * t243;
t244 = cos(t248);
t276 = pkin(2) * t244;
t250 = cos(pkin(9));
t275 = t250 * pkin(4);
t274 = pkin(1) * qJD(1);
t246 = pkin(9) + qJ(5);
t240 = sin(t246);
t273 = Icges(6,4) * t240;
t241 = cos(t246);
t272 = Icges(6,4) * t241;
t242 = pkin(8) + t248;
t237 = sin(t242);
t270 = qJD(5) * t237;
t238 = cos(t242);
t269 = qJD(5) * t238;
t252 = sin(qJ(1));
t268 = t252 * t274;
t253 = cos(qJ(1));
t267 = t253 * t274;
t266 = -t238 * pkin(3) - t237 * qJ(4) - t276;
t265 = qJD(4) * t238 - t267;
t249 = sin(pkin(9));
t264 = -rSges(5,1) * t250 + rSges(5,2) * t249;
t263 = rSges(6,1) * t241 - rSges(6,2) * t240;
t262 = Icges(6,1) * t241 - t273;
t261 = -Icges(6,2) * t240 + t272;
t260 = Icges(6,5) * t241 - Icges(6,6) * t240;
t219 = Icges(6,6) * t238 - t261 * t237;
t221 = Icges(6,5) * t238 - t262 * t237;
t259 = -t219 * t240 + t221 * t241;
t220 = Icges(6,6) * t237 + t261 * t238;
t222 = Icges(6,5) * t237 + t262 * t238;
t258 = t220 * t240 - t222 * t241;
t230 = Icges(6,2) * t241 + t273;
t231 = Icges(6,1) * t240 + t272;
t257 = t230 * t240 - t231 * t241;
t247 = qJD(1) + qJD(2);
t256 = t247 * (-t237 * pkin(3) + t238 * qJ(4)) + qJD(4) * t237 - t268;
t254 = qJD(3) ^ 2;
t234 = t253 * rSges(2,1) - t252 * rSges(2,2);
t233 = -t252 * rSges(2,1) - t253 * rSges(2,2);
t232 = t240 * rSges(6,1) + t241 * rSges(6,2);
t229 = Icges(6,5) * t240 + Icges(6,6) * t241;
t226 = -t267 - t247 * (t244 * rSges(3,1) - t243 * rSges(3,2));
t225 = -t268 + t247 * (-t243 * rSges(3,1) - t244 * rSges(3,2));
t224 = t237 * rSges(6,3) + t263 * t238;
t223 = t238 * rSges(6,3) - t263 * t237;
t218 = Icges(6,3) * t237 + t260 * t238;
t217 = Icges(6,3) * t238 - t260 * t237;
t216 = -t267 + (-t238 * rSges(4,1) + t237 * rSges(4,2) - t276) * t247;
t215 = -t268 + (-t237 * rSges(4,1) - t238 * rSges(4,2) - t277) * t247;
t214 = (-t237 * rSges(5,3) + t264 * t238 + t266) * t247 + t265;
t213 = (t238 * rSges(5,3) + t264 * t237 - t277) * t247 + t256;
t212 = qJD(3) + (-t223 * t237 + t224 * t238) * qJD(5);
t211 = t232 * t270 + (-pkin(7) * t237 - t275 * t238 - t224 + t266) * t247 + t265;
t210 = -t232 * t269 + (pkin(7) * t238 - t275 * t237 + t223 - t277) * t247 + t256;
t1 = m(3) * (t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(4) * (t215 ^ 2 + t216 ^ 2 + t254) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t254) / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t247 * ((t241 * t230 + t240 * t231) * t247 + ((t241 * t219 + t240 * t221) * t238 + (t241 * t220 + t240 * t222) * t237) * qJD(5)) / 0.2e1 + ((t238 * t229 + t257 * t237) * t247 + (t238 ^ 2 * t217 + (t258 * t237 + (t218 - t259) * t238) * t237) * qJD(5)) * t269 / 0.2e1 + ((t237 * t229 - t257 * t238) * t247 + (t237 ^ 2 * t218 + (t259 * t238 + (t217 - t258) * t237) * t238) * qJD(5)) * t270 / 0.2e1 + (m(2) * (t233 ^ 2 + t234 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t250 ^ 2 + (Icges(5,1) * t249 + 0.2e1 * Icges(5,4) * t250) * t249) * t247 ^ 2 / 0.2e1;
T = t1;
