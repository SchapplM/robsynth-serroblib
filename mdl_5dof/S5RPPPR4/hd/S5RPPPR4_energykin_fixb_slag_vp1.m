% Calculate kinetic energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:05
% EndTime: 2019-12-31 17:45:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (364->99), mult. (325->152), div. (0->0), fcn. (238->8), ass. (0->56)
t251 = qJ(1) + pkin(7);
t247 = sin(t251);
t249 = cos(t251);
t277 = qJD(1) * t249 * qJ(4) + qJD(4) * t247;
t252 = sin(pkin(8));
t275 = pkin(4) * t252;
t255 = sin(qJ(1));
t274 = t255 * pkin(1);
t250 = pkin(8) + qJ(5);
t246 = sin(t250);
t273 = Icges(6,4) * t246;
t248 = cos(t250);
t272 = Icges(6,4) * t248;
t256 = cos(qJ(1));
t245 = qJD(1) * t256 * pkin(1);
t271 = qJD(1) * (t249 * pkin(2) + t247 * qJ(3)) + t245;
t244 = qJD(3) * t247;
t270 = qJD(4) * t249 + t244;
t269 = qJD(5) * t247;
t268 = -t247 * pkin(2) + t249 * qJ(3) - t274;
t267 = -qJD(3) * t249 + t271;
t253 = cos(pkin(8));
t266 = rSges(5,1) * t252 + rSges(5,2) * t253;
t265 = rSges(6,1) * t246 + rSges(6,2) * t248;
t264 = Icges(6,1) * t246 + t272;
t263 = Icges(6,2) * t248 + t273;
t262 = Icges(6,5) * t246 + Icges(6,6) * t248;
t225 = Icges(6,6) * t249 + t263 * t247;
t227 = Icges(6,5) * t249 + t264 * t247;
t261 = -t225 * t248 - t227 * t246;
t226 = Icges(6,6) * t247 - t263 * t249;
t228 = Icges(6,5) * t247 - t264 * t249;
t260 = t226 * t248 + t228 * t246;
t235 = -Icges(6,2) * t246 + t272;
t236 = Icges(6,1) * t248 - t273;
t259 = t235 * t248 + t236 * t246;
t257 = qJD(2) ^ 2;
t254 = -pkin(6) - qJ(4);
t240 = t256 * rSges(2,1) - t255 * rSges(2,2);
t239 = t255 * rSges(2,1) + t256 * rSges(2,2);
t238 = t248 * rSges(6,1) - t246 * rSges(6,2);
t234 = Icges(6,5) * t248 - Icges(6,6) * t246;
t232 = t245 + qJD(1) * (t249 * rSges(3,1) - t247 * rSges(3,2));
t231 = (-t247 * rSges(3,1) - t249 * rSges(3,2) - t274) * qJD(1);
t230 = t247 * rSges(6,3) - t265 * t249;
t229 = t249 * rSges(6,3) + t265 * t247;
t224 = Icges(6,3) * t247 - t262 * t249;
t223 = Icges(6,3) * t249 + t262 * t247;
t222 = qJD(1) * (-t249 * rSges(4,2) + t247 * rSges(4,3)) + t267;
t221 = t244 + (t247 * rSges(4,2) + t249 * rSges(4,3) + t268) * qJD(1);
t220 = qJD(1) * (t249 * rSges(5,3) + t266 * t247) + t267 + t277;
t219 = (t266 * t249 + (-rSges(5,3) - qJ(4)) * t247 + t268) * qJD(1) + t270;
t218 = qJD(2) + (-t229 * t247 + t230 * t249) * qJD(5);
t217 = (t247 * t275 + t229) * qJD(1) + (-qJD(3) + qJD(1) * (-qJ(4) - t254) - qJD(5) * t238) * t249 + t271 + t277;
t216 = t238 * t269 + (t247 * t254 + t249 * t275 - t230 + t268) * qJD(1) + t270;
t1 = m(3) * (t231 ^ 2 + t232 ^ 2 + t257) / 0.2e1 + m(4) * (t221 ^ 2 + t222 ^ 2 + t257) / 0.2e1 + m(5) * (t219 ^ 2 + t220 ^ 2 + t257) / 0.2e1 + m(6) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + qJD(5) * t249 * ((t249 * t234 + t259 * t247) * qJD(1) + (t249 ^ 2 * t223 + (t260 * t247 + (t224 - t261) * t249) * t247) * qJD(5)) / 0.2e1 + ((t247 * t234 - t259 * t249) * qJD(1) + (t247 ^ 2 * t224 + (t261 * t249 + (t223 - t260) * t247) * t249) * qJD(5)) * t269 / 0.2e1 + qJD(1) * ((-t246 * t235 + t248 * t236) * qJD(1) + ((-t246 * t225 + t248 * t227) * t249 + (-t246 * t226 + t248 * t228) * t247) * qJD(5)) / 0.2e1 + (m(2) * (t239 ^ 2 + t240 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1) + Icges(5,1) * t253 ^ 2 + (-0.2e1 * Icges(5,4) * t253 + Icges(5,2) * t252) * t252) * qJD(1) ^ 2 / 0.2e1;
T = t1;
