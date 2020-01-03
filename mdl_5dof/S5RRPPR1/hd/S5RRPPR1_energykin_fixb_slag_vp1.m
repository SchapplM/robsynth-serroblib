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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:55:42
% EndTime: 2020-01-03 11:55:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (500->94), mult. (318->156), div. (0->0), fcn. (232->10), ass. (0->59)
t249 = qJD(1) + qJD(2);
t251 = sin(pkin(9));
t252 = cos(pkin(9));
t278 = -qJD(4) + t249 * (rSges(5,1) * t252 - rSges(5,2) * t251);
t275 = pkin(2) * t249;
t274 = t252 * pkin(4);
t273 = pkin(1) * qJD(1);
t248 = pkin(9) + qJ(5);
t242 = sin(t248);
t272 = Icges(6,4) * t242;
t243 = cos(t248);
t271 = Icges(6,4) * t243;
t250 = qJ(1) + qJ(2);
t254 = sin(qJ(1));
t240 = t254 * t273;
t245 = sin(t250);
t269 = t245 * t275 + t240;
t255 = cos(qJ(1));
t241 = t255 * t273;
t246 = cos(t250);
t268 = t246 * t275 + t241;
t244 = pkin(8) + t250;
t237 = sin(t244);
t267 = qJD(5) * t237;
t238 = cos(t244);
t266 = t249 * (t237 * pkin(3) - t238 * qJ(4)) + t269;
t264 = rSges(6,1) * t243 - rSges(6,2) * t242;
t263 = Icges(6,1) * t243 - t272;
t262 = -Icges(6,2) * t242 + t271;
t261 = Icges(6,5) * t243 - Icges(6,6) * t242;
t219 = -Icges(6,6) * t238 + t237 * t262;
t221 = -Icges(6,5) * t238 + t237 * t263;
t260 = -t219 * t242 + t221 * t243;
t220 = -Icges(6,6) * t237 - t238 * t262;
t222 = -Icges(6,5) * t237 - t238 * t263;
t259 = t220 * t242 - t222 * t243;
t230 = Icges(6,2) * t243 + t272;
t231 = Icges(6,1) * t242 + t271;
t258 = t230 * t242 - t231 * t243;
t256 = qJD(3) ^ 2;
t236 = -t255 * rSges(2,1) + t254 * rSges(2,2);
t235 = t254 * rSges(2,1) + t255 * rSges(2,2);
t232 = t242 * rSges(6,1) + t243 * rSges(6,2);
t229 = Icges(6,5) * t242 + Icges(6,6) * t243;
t228 = -t238 * pkin(3) - t237 * qJ(4);
t226 = t241 - t249 * (-t246 * rSges(3,1) + t245 * rSges(3,2));
t225 = t240 + t249 * (t245 * rSges(3,1) + t246 * rSges(3,2));
t224 = -t237 * rSges(6,3) - t238 * t264;
t223 = -t238 * rSges(6,3) + t237 * t264;
t218 = -Icges(6,3) * t237 - t238 * t261;
t217 = -Icges(6,3) * t238 + t237 * t261;
t216 = -t249 * (-t238 * rSges(4,1) + t237 * rSges(4,2)) + t268;
t215 = t249 * (t237 * rSges(4,1) + t238 * rSges(4,2)) + t269;
t214 = (t237 * rSges(5,3) - t228) * t249 + t278 * t238 + t268;
t213 = -t249 * t238 * rSges(5,3) + t278 * t237 + t266;
t212 = qJD(3) + (t223 * t237 - t224 * t238) * qJD(5);
t211 = -t232 * t267 - qJD(4) * t238 + (pkin(7) * t237 + t238 * t274 - t224 - t228) * t249 + t268;
t210 = t249 * t223 + (-pkin(7) * t249 + qJD(5) * t232) * t238 + (t249 * t274 - qJD(4)) * t237 + t266;
t1 = m(3) * (t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(4) * (t215 ^ 2 + t216 ^ 2 + t256) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t256) / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t249 * ((t243 * t230 + t242 * t231) * t249 + (-(t243 * t219 + t242 * t221) * t238 - (t243 * t220 + t242 * t222) * t237) * qJD(5)) / 0.2e1 - qJD(5) * t238 * ((-t238 * t229 - t237 * t258) * t249 + (t238 ^ 2 * t217 + (t259 * t237 + (t218 - t260) * t238) * t237) * qJD(5)) / 0.2e1 - ((-t237 * t229 + t238 * t258) * t249 + (t237 ^ 2 * t218 + (t260 * t238 + (t217 - t259) * t237) * t238) * qJD(5)) * t267 / 0.2e1 + (m(2) * (t235 ^ 2 + t236 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t252 ^ 2 + (Icges(5,1) * t251 + 0.2e1 * Icges(5,4) * t252) * t251) * t249 ^ 2 / 0.2e1;
T = t1;
