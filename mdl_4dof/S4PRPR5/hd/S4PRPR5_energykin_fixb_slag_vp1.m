% Calculate kinetic energy for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:57
% EndTime: 2019-12-31 16:22:57
% DurationCPUTime: 0.72s
% Computational Cost: add. (494->102), mult. (731->179), div. (0->0), fcn. (708->8), ass. (0->63)
t298 = Icges(3,3) + Icges(4,3);
t246 = qJ(2) + pkin(7);
t244 = sin(t246);
t245 = cos(t246);
t251 = sin(qJ(2));
t253 = cos(qJ(2));
t297 = Icges(3,5) * t253 + Icges(4,5) * t245 - Icges(3,6) * t251 - Icges(4,6) * t244;
t247 = sin(pkin(6));
t248 = cos(pkin(6));
t296 = t247 * t248;
t295 = t297 * t247 - t298 * t248;
t294 = t298 * t247 + t297 * t248;
t288 = t248 ^ 2;
t289 = t247 ^ 2;
t291 = t288 + t289;
t290 = t291 * qJD(2);
t287 = qJD(2) ^ 2;
t286 = pkin(2) * t251;
t285 = t253 * pkin(2);
t283 = t244 * t247;
t282 = t244 * t248;
t250 = sin(qJ(4));
t281 = t247 * t250;
t252 = cos(qJ(4));
t280 = t247 * t252;
t279 = t248 * t250;
t278 = t248 * t252;
t277 = qJD(2) * t247;
t276 = qJD(2) * t248;
t275 = qJD(3) * t248;
t274 = qJD(4) * t244;
t273 = qJD(4) * t245;
t272 = qJD(1) + (-qJ(3) * t248 + t285 * t247) * t277 + (qJ(3) * t247 + t285 * t248) * t276;
t268 = qJD(2) * (-t244 * rSges(4,1) - t245 * rSges(4,2) - t286);
t267 = qJD(2) * (-t244 * pkin(3) + t245 * pkin(5) - t286);
t243 = qJD(3) * t247;
t241 = t251 * rSges(3,1) + t253 * rSges(3,2);
t237 = t247 * t274 - t276;
t236 = t248 * t274 + t277;
t235 = t245 * t278 + t281;
t234 = -t245 * t279 + t280;
t233 = t245 * t280 - t279;
t232 = -t245 * t281 - t278;
t219 = -t245 * rSges(5,3) + (rSges(5,1) * t252 - rSges(5,2) * t250) * t244;
t218 = -Icges(5,5) * t245 + (Icges(5,1) * t252 - Icges(5,4) * t250) * t244;
t217 = -Icges(5,6) * t245 + (Icges(5,4) * t252 - Icges(5,2) * t250) * t244;
t216 = -Icges(5,3) * t245 + (Icges(5,5) * t252 - Icges(5,6) * t250) * t244;
t213 = t248 * t268 + t243;
t212 = t247 * t268 - t275;
t211 = t235 * rSges(5,1) + t234 * rSges(5,2) + rSges(5,3) * t282;
t210 = t233 * rSges(5,1) + t232 * rSges(5,2) + rSges(5,3) * t283;
t209 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t282;
t208 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t283;
t207 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t282;
t206 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t283;
t205 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t282;
t204 = Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t283;
t203 = qJD(1) + (rSges(3,1) * t253 - rSges(3,2) * t251) * t290;
t202 = t272 + (rSges(4,1) * t245 - rSges(4,2) * t244) * t290;
t201 = t210 * t273 + t237 * t219 + t248 * t267 + t243;
t200 = -t211 * t273 - t236 * t219 + t247 * t267 - t275;
t199 = t236 * t210 - t237 * t211 + t272 + (pkin(3) * t245 + pkin(5) * t244) * t290;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t291 * t287 * t241 ^ 2 + t203 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(5) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + t236 * ((t205 * t282 + t234 * t207 + t235 * t209) * t236 + (t204 * t282 + t234 * t206 + t235 * t208) * t237 - (t216 * t282 + t234 * t217 + t235 * t218) * t273) / 0.2e1 + t237 * ((t205 * t283 + t232 * t207 + t233 * t209) * t236 + (t204 * t283 + t232 * t206 + t233 * t208) * t237 - (t216 * t283 + t232 * t217 + t233 * t218) * t273) / 0.2e1 - ((-t204 * t237 - t205 * t236 + t216 * t273) * t245 + ((-t207 * t250 + t209 * t252) * t236 + (-t206 * t250 + t208 * t252) * t237 - (-t217 * t250 + t218 * t252) * t273) * t244) * t273 / 0.2e1 + (t294 * t289 - t295 * t296) * t247 * t287 / 0.2e1 - (t295 * t288 - t294 * t296) * t248 * t287 / 0.2e1;
T = t1;
