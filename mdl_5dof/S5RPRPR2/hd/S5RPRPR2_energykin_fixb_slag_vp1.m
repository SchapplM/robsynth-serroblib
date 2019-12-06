% Calculate kinetic energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:19
% EndTime: 2019-12-05 17:49:19
% DurationCPUTime: 0.41s
% Computational Cost: add. (492->94), mult. (319->154), div. (0->0), fcn. (232->10), ass. (0->59)
t248 = qJ(1) + pkin(8);
t244 = qJ(3) + t248;
t238 = cos(t244);
t277 = qJD(5) * t238;
t252 = sin(qJ(1));
t274 = pkin(1) * t252;
t253 = cos(qJ(1));
t273 = pkin(1) * t253;
t250 = cos(pkin(9));
t272 = pkin(4) * t250;
t246 = pkin(9) + qJ(5);
t240 = sin(t246);
t271 = Icges(6,4) * t240;
t242 = cos(t246);
t270 = Icges(6,4) * t242;
t237 = sin(t244);
t268 = qJD(5) * t237;
t249 = sin(pkin(9));
t267 = -rSges(5,1) * t250 + rSges(5,2) * t249;
t266 = rSges(6,1) * t242 - rSges(6,2) * t240;
t265 = Icges(6,1) * t242 - t271;
t264 = -Icges(6,2) * t240 + t270;
t263 = Icges(6,5) * t242 - Icges(6,6) * t240;
t219 = Icges(6,6) * t238 - t264 * t237;
t221 = Icges(6,5) * t238 - t265 * t237;
t262 = -t219 * t240 + t221 * t242;
t220 = Icges(6,6) * t237 + t264 * t238;
t222 = Icges(6,5) * t237 + t265 * t238;
t261 = t220 * t240 - t222 * t242;
t230 = Icges(6,2) * t242 + t271;
t231 = Icges(6,1) * t240 + t270;
t260 = t230 * t240 - t231 * t242;
t241 = sin(t248);
t259 = (-pkin(2) * t241 - t274) * qJD(1);
t243 = cos(t248);
t258 = (-pkin(2) * t243 - t273) * qJD(1);
t257 = qJD(4) * t238 + t258;
t247 = qJD(1) + qJD(3);
t256 = t247 * (-pkin(3) * t237 + qJ(4) * t238) + qJD(4) * t237 + t259;
t254 = qJD(2) ^ 2;
t234 = rSges(2,1) * t253 - rSges(2,2) * t252;
t233 = -rSges(2,1) * t252 - rSges(2,2) * t253;
t232 = rSges(6,1) * t240 + rSges(6,2) * t242;
t229 = Icges(6,5) * t240 + Icges(6,6) * t242;
t228 = t238 * pkin(3) + t237 * qJ(4);
t226 = (-rSges(3,1) * t243 + rSges(3,2) * t241 - t273) * qJD(1);
t225 = (-rSges(3,1) * t241 - rSges(3,2) * t243 - t274) * qJD(1);
t224 = rSges(6,3) * t237 + t266 * t238;
t223 = rSges(6,3) * t238 - t266 * t237;
t218 = Icges(6,3) * t237 + t263 * t238;
t217 = Icges(6,3) * t238 - t263 * t237;
t216 = -t247 * (rSges(4,1) * t238 - rSges(4,2) * t237) + t258;
t215 = t247 * (-rSges(4,1) * t237 - rSges(4,2) * t238) + t259;
t214 = (-t237 * rSges(5,3) + t267 * t238 - t228) * t247 + t257;
t213 = t247 * (rSges(5,3) * t238 + t267 * t237) + t256;
t212 = qJD(2) + (-t223 * t237 + t224 * t238) * qJD(5);
t211 = t232 * t268 + (-pkin(7) * t237 - t272 * t238 - t224 - t228) * t247 + t257;
t210 = -t232 * t277 + t256 + (pkin(7) * t238 - t272 * t237 + t223) * t247;
t1 = m(3) * (t225 ^ 2 + t226 ^ 2 + t254) / 0.2e1 + m(4) * (t215 ^ 2 + t216 ^ 2 + t254) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t254) / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t247 * ((t242 * t230 + t240 * t231) * t247 + ((t219 * t242 + t221 * t240) * t238 + (t220 * t242 + t222 * t240) * t237) * qJD(5)) / 0.2e1 + ((t238 * t229 + t260 * t237) * t247 + (t238 ^ 2 * t217 + (t261 * t237 + (t218 - t262) * t238) * t237) * qJD(5)) * t277 / 0.2e1 + ((t237 * t229 - t260 * t238) * t247 + (t237 ^ 2 * t218 + (t262 * t238 + (t217 - t261) * t237) * t238) * qJD(5)) * t268 / 0.2e1 + (Icges(4,3) + Icges(5,2) * t250 ^ 2 + (Icges(5,1) * t249 + 0.2e1 * Icges(5,4) * t250) * t249) * t247 ^ 2 / 0.2e1 + (m(2) * (t233 ^ 2 + t234 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
