% Calculate kinetic energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:41
% EndTime: 2019-12-31 19:05:42
% DurationCPUTime: 0.82s
% Computational Cost: add. (573->111), mult. (977->193), div. (0->0), fcn. (1096->8), ass. (0->68)
t309 = sin(qJ(3));
t310 = sin(qJ(1));
t311 = cos(qJ(3));
t312 = cos(qJ(1));
t259 = -t310 * t309 - t312 * t311;
t319 = t259 ^ 2;
t260 = t312 * t309 - t310 * t311;
t318 = t260 ^ 2;
t279 = sin(qJ(4));
t280 = cos(qJ(4));
t317 = -Icges(5,5) * t279 - Icges(5,6) * t280;
t316 = t260 * t259;
t277 = qJD(1) - qJD(3);
t315 = t317 * t277;
t308 = t280 * pkin(4);
t278 = qJ(4) + qJ(5);
t275 = sin(t278);
t304 = Icges(6,4) * t275;
t276 = cos(t278);
t303 = Icges(6,4) * t276;
t302 = qJD(4) * (-t279 * rSges(5,1) - t280 * rSges(5,2));
t301 = qJD(4) + qJD(5);
t300 = pkin(4) * qJD(4) * t279;
t299 = -qJD(2) * t312 + qJD(1) * (t312 * pkin(1) + t310 * qJ(2));
t298 = -rSges(5,1) * t280 + rSges(5,2) * t279;
t297 = -rSges(6,1) * t276 + rSges(6,2) * t275;
t295 = -Icges(6,1) * t276 + t304;
t293 = Icges(6,2) * t275 - t303;
t292 = -Icges(5,5) * t280 + Icges(5,6) * t279;
t291 = -Icges(6,5) * t276 + Icges(6,6) * t275;
t287 = qJD(1) * t312 * pkin(2) + t299;
t286 = t277 * (-t259 * pkin(3) + t260 * pkin(7)) + t287;
t266 = t310 * pkin(1) - t312 * qJ(2);
t274 = qJD(2) * t310;
t285 = t274 + (-t310 * pkin(2) - t266) * qJD(1);
t250 = t301 * t260;
t251 = t301 * t259;
t284 = -(-Icges(6,3) * t259 + t291 * t260) * t251 + (Icges(6,3) * t260 + t291 * t259) * t250 + (-Icges(6,5) * t275 - Icges(6,6) * t276) * t277;
t235 = -Icges(6,6) * t259 + t293 * t260;
t236 = Icges(6,6) * t260 + t293 * t259;
t237 = -Icges(6,5) * t259 + t295 * t260;
t238 = Icges(6,5) * t260 + t295 * t259;
t256 = -Icges(6,2) * t276 - t304;
t257 = -Icges(6,1) * t275 - t303;
t283 = (t236 * t275 - t238 * t276) * t250 - (t235 * t275 - t237 * t276) * t251 + (t256 * t275 - t257 * t276) * t277;
t268 = t312 * rSges(2,1) - t310 * rSges(2,2);
t267 = t310 * rSges(2,1) + t312 * rSges(2,2);
t258 = -t275 * rSges(6,1) - t276 * rSges(6,2);
t254 = qJD(1) * (t312 * rSges(3,1) + t310 * rSges(3,3)) + t299;
t253 = t274 + (-t310 * rSges(3,1) + t312 * rSges(3,3) - t266) * qJD(1);
t252 = -t260 * pkin(3) - t259 * pkin(7);
t248 = t260 * rSges(5,3) + t298 * t259;
t247 = -t259 * rSges(5,3) + t298 * t260;
t242 = Icges(5,3) * t260 + t292 * t259;
t241 = -Icges(5,3) * t259 + t292 * t260;
t240 = t260 * rSges(6,3) + t297 * t259;
t239 = -t259 * rSges(6,3) + t297 * t260;
t232 = t277 * (-t259 * rSges(4,1) - t260 * rSges(4,2)) + t287;
t231 = -t277 * (-t260 * rSges(4,1) + t259 * rSges(4,2)) + t285;
t230 = pkin(8) * t260 - t308 * t259;
t229 = -pkin(8) * t259 - t308 * t260;
t228 = (t247 * t260 + t248 * t259) * qJD(4);
t227 = t277 * t248 - t260 * t302 + t286;
t226 = -t259 * t302 + (-t247 - t252) * t277 + t285;
t225 = t260 * t300 - t250 * t258 + (t230 + t240) * t277 + t286;
t224 = t259 * t300 - t251 * t258 + (-t229 - t239 - t252) * t277 + t285;
t223 = t250 * t239 + t251 * t240 + (t229 * t260 + t230 * t259) * qJD(4);
t1 = m(3) * (t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (t231 ^ 2 + t232 ^ 2) / 0.2e1 + t277 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + qJD(4) * t260 * (t260 * t315 + (-t241 * t316 + t318 * t242) * qJD(4)) / 0.2e1 - qJD(4) * t259 * (-t259 * t315 + (t319 * t241 - t242 * t316) * qJD(4)) / 0.2e1 + m(6) * (t223 ^ 2 + t224 ^ 2 + t225 ^ 2) / 0.2e1 + t250 * (t283 * t259 + t284 * t260) / 0.2e1 - t251 * (-t284 * t259 + t283 * t260) / 0.2e1 + ((t318 + t319) * qJD(4) * t317 + (-t276 * t236 - t275 * t238) * t250 - (-t276 * t235 - t275 * t237) * t251 + (t280 ^ 2 * Icges(5,2) - t276 * t256 - t275 * t257 + (Icges(5,1) * t279 + 0.2e1 * Icges(5,4) * t280) * t279) * t277) * t277 / 0.2e1 + (m(2) * (t267 ^ 2 + t268 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
