% Calculate kinetic energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:03
% EndTime: 2019-12-31 18:01:04
% DurationCPUTime: 0.48s
% Computational Cost: add. (461->80), mult. (528->131), div. (0->0), fcn. (544->8), ass. (0->51)
t280 = pkin(8) + qJ(4);
t274 = sin(t280);
t275 = cos(t280);
t284 = sin(qJ(1));
t285 = cos(qJ(1));
t237 = -t274 * t284 - t275 * t285;
t292 = t237 ^ 2;
t238 = t274 * t285 - t275 * t284;
t291 = t238 ^ 2;
t258 = sin(qJ(5));
t259 = cos(qJ(5));
t290 = -Icges(6,5) * t258 - Icges(6,6) * t259;
t288 = t238 * t237;
t255 = qJD(1) - qJD(4);
t287 = t290 * t255;
t281 = qJD(5) * (-t258 * rSges(6,1) - t259 * rSges(6,2));
t279 = t285 * pkin(2);
t278 = t284 * pkin(2);
t256 = sin(pkin(8));
t277 = t285 * t256;
t276 = t284 * t256;
t246 = pkin(1) * t284 - qJ(2) * t285;
t273 = -t246 - t278;
t272 = -qJD(2) * t285 + qJD(1) * (pkin(1) * t285 + qJ(2) * t284);
t271 = -rSges(6,1) * t259 + rSges(6,2) * t258;
t268 = -Icges(6,5) * t259 + Icges(6,6) * t258;
t264 = qJD(1) * t279 + t272;
t257 = cos(pkin(8));
t252 = t257 * pkin(3) + pkin(2);
t263 = qJD(1) * (pkin(3) * t276 + t252 * t285 - t279) + t264;
t254 = qJD(2) * t284;
t262 = t254 + (pkin(3) * t277 - t252 * t284 + t273 + t278) * qJD(1);
t260 = qJD(3) ^ 2;
t248 = rSges(2,1) * t285 - rSges(2,2) * t284;
t247 = rSges(2,1) * t284 + rSges(2,2) * t285;
t240 = t257 * t284 - t277;
t239 = -t257 * t285 - t276;
t234 = qJD(1) * (rSges(3,1) * t285 + rSges(3,3) * t284) + t272;
t233 = t254 + (-rSges(3,1) * t284 + rSges(3,3) * t285 - t246) * qJD(1);
t232 = qJD(1) * (-t239 * rSges(4,1) + t240 * rSges(4,2)) + t264;
t231 = t254 + (-t240 * rSges(4,1) - t239 * rSges(4,2) + t273) * qJD(1);
t230 = t238 * rSges(6,3) + t237 * t271;
t229 = -t237 * rSges(6,3) + t238 * t271;
t224 = Icges(6,3) * t238 + t237 * t268;
t223 = -Icges(6,3) * t237 + t238 * t268;
t222 = t255 * (-t237 * rSges(5,1) - t238 * rSges(5,2)) + t263;
t221 = -t255 * (-t238 * rSges(5,1) + t237 * rSges(5,2)) + t262;
t220 = -qJD(3) + (t229 * t238 + t230 * t237) * qJD(5);
t219 = -t238 * t281 + (-t237 * pkin(4) + t238 * pkin(7) + t230) * t255 + t263;
t218 = -t237 * t281 + (t238 * pkin(4) + t237 * pkin(7) - t229) * t255 + t262;
t1 = m(3) * (t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(4) * (t231 ^ 2 + t232 ^ 2 + t260) / 0.2e1 + m(5) * (t221 ^ 2 + t222 ^ 2 + t260) / 0.2e1 + t255 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + qJD(5) * t238 * (t238 * t287 + (-t223 * t288 + t291 * t224) * qJD(5)) / 0.2e1 - qJD(5) * t237 * (-t237 * t287 + (t292 * t223 - t224 * t288) * qJD(5)) / 0.2e1 + t255 * ((t259 ^ 2 * Icges(6,2) + (Icges(6,1) * t258 + 0.2e1 * Icges(6,4) * t259) * t258) * t255 + (t291 + t292) * t290 * qJD(5)) / 0.2e1 + (m(2) * (t247 ^ 2 + t248 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
