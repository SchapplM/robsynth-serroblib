% Calculate kinetic energy for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:17
% EndTime: 2019-12-05 16:17:18
% DurationCPUTime: 0.32s
% Computational Cost: add. (591->101), mult. (503->172), div. (0->0), fcn. (486->8), ass. (0->52)
t271 = pkin(2) * qJD(2);
t250 = pkin(8) + qJ(2);
t248 = qJ(3) + t250;
t244 = sin(t248);
t252 = sin(pkin(9));
t270 = t244 * t252;
t245 = cos(t248);
t269 = t245 * t252;
t253 = cos(pkin(9));
t254 = sin(qJ(5));
t268 = t253 * t254;
t255 = cos(qJ(5));
t267 = t253 * t255;
t247 = cos(t250);
t243 = t247 * t271;
t251 = qJD(2) + qJD(3);
t266 = t251 * (pkin(3) * t245 + qJ(4) * t244) + t243;
t265 = qJD(5) * t252;
t246 = sin(t250);
t264 = t246 * t271;
t262 = qJD(4) * t244 - t264;
t261 = rSges(5,1) * t253 - rSges(5,2) * t252;
t229 = -t244 * t268 - t245 * t255;
t230 = t244 * t267 - t245 * t254;
t231 = t244 * t255 - t245 * t268;
t232 = t244 * t254 + t245 * t267;
t260 = (Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t270) * t244 + (Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t269) * t245;
t259 = t260 * t252;
t258 = -t251 * (pkin(4) * t253 + pkin(7) * t252) + (-rSges(6,3) * t253 + (rSges(6,1) * t255 - rSges(6,2) * t254) * t252) * t265;
t257 = qJD(1) ^ 2;
t256 = qJD(2) ^ 2;
t241 = -t253 * qJD(5) + t251;
t240 = rSges(3,1) * t247 - rSges(3,2) * t246;
t239 = rSges(3,1) * t246 + rSges(3,2) * t247;
t238 = pkin(3) * t244 - qJ(4) * t245;
t235 = -Icges(6,5) * t253 + (Icges(6,1) * t255 - Icges(6,4) * t254) * t252;
t234 = -Icges(6,6) * t253 + (Icges(6,4) * t255 - Icges(6,2) * t254) * t252;
t233 = -Icges(6,3) * t253 + (Icges(6,5) * t255 - Icges(6,6) * t254) * t252;
t228 = t243 + t251 * (rSges(4,1) * t245 - rSges(4,2) * t244);
t227 = -t264 - t251 * (rSges(4,1) * t244 + rSges(4,2) * t245);
t226 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t269;
t225 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t270;
t224 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t269;
t223 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t270;
t222 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t269;
t221 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t270;
t218 = t251 * t244 * rSges(5,3) + (t251 * t261 - qJD(4)) * t245 + t266;
t217 = (t245 * rSges(5,3) - t261 * t244 - t238) * t251 + t262;
t216 = qJD(1) + (t225 * t245 - t226 * t244) * t265;
t215 = t241 * t226 + (-qJD(4) - t258) * t245 + t266;
t214 = -t241 * t225 - t251 * t238 + t258 * t244 + t262;
t1 = m(2) * t257 / 0.2e1 + m(3) * (t257 + (t239 ^ 2 + t240 ^ 2) * t256) / 0.2e1 + t256 * Icges(3,3) / 0.2e1 + m(4) * (t227 ^ 2 + t228 ^ 2 + t257) / 0.2e1 + m(5) * (t217 ^ 2 + t218 ^ 2 + t257) / 0.2e1 + m(6) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + t241 * ((-t253 * t233 + (-t234 * t254 + t235 * t255) * t252) * t241 + (((-t222 * t254 + t224 * t255) * t245 + (-t221 * t254 + t223 * t255) * t244) * t252 - t260 * t253) * t265) / 0.2e1 + (Icges(4,3) + Icges(5,2) * t253 ^ 2 + (Icges(5,1) * t252 + 0.2e1 * Icges(5,4) * t253) * t252) * t251 ^ 2 / 0.2e1 + (t245 * ((t231 * t234 + t232 * t235 + t233 * t269) * t241 + ((t221 * t231 + t223 * t232) * t244 + (t231 * t222 + t232 * t224 + t259) * t245) * t265) + t244 * ((t229 * t234 + t230 * t235 + t233 * t270) * t241 + ((t222 * t229 + t224 * t230) * t245 + (t229 * t221 + t230 * t223 + t259) * t244) * t265)) * t265 / 0.2e1;
T = t1;
