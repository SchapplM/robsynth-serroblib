% Calculate kinetic energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:57
% EndTime: 2022-01-20 10:33:57
% DurationCPUTime: 0.30s
% Computational Cost: add. (482->87), mult. (292->145), div. (0->0), fcn. (206->10), ass. (0->58)
t232 = qJ(1) + qJ(2);
t228 = sin(t232);
t255 = pkin(2) * t228;
t254 = pkin(1) * qJD(1);
t233 = sin(qJ(5));
t253 = Icges(6,4) * t233;
t235 = cos(qJ(5));
t252 = Icges(6,4) * t235;
t236 = cos(qJ(1));
t225 = t236 * t254;
t229 = cos(t232);
t231 = qJD(1) + qJD(2);
t251 = t231 * pkin(2) * t229 + t225;
t227 = pkin(9) + t232;
t224 = qJ(4) + t227;
t220 = sin(t224);
t250 = qJD(5) * t220;
t221 = cos(t224);
t249 = qJD(5) * t221;
t234 = sin(qJ(1));
t248 = t234 * t254;
t223 = cos(t227);
t247 = t231 * pkin(3) * t223 + t251;
t246 = rSges(6,1) * t235 - rSges(6,2) * t233;
t245 = Icges(6,1) * t235 - t253;
t244 = -Icges(6,2) * t233 + t252;
t243 = Icges(6,5) * t235 - Icges(6,6) * t233;
t204 = -Icges(6,6) * t221 + t244 * t220;
t206 = -Icges(6,5) * t221 + t245 * t220;
t242 = t204 * t233 - t206 * t235;
t205 = Icges(6,6) * t220 + t244 * t221;
t207 = Icges(6,5) * t220 + t245 * t221;
t241 = -t205 * t233 + t207 * t235;
t214 = Icges(6,2) * t235 + t253;
t215 = Icges(6,1) * t233 + t252;
t240 = -t214 * t233 + t215 * t235;
t222 = sin(t227);
t239 = (-pkin(3) * t222 - t255) * t231 - t248;
t237 = qJD(3) ^ 2;
t226 = qJD(4) + t231;
t219 = rSges(2,1) * t236 - rSges(2,2) * t234;
t218 = rSges(2,1) * t234 + rSges(2,2) * t236;
t217 = rSges(6,1) * t233 + rSges(6,2) * t235;
t213 = Icges(6,5) * t233 + Icges(6,6) * t235;
t211 = t225 + t231 * (rSges(3,1) * t229 - rSges(3,2) * t228);
t210 = -t248 - t231 * (rSges(3,1) * t228 + rSges(3,2) * t229);
t209 = rSges(6,3) * t220 + t246 * t221;
t208 = -rSges(6,3) * t221 + t246 * t220;
t203 = Icges(6,3) * t220 + t243 * t221;
t202 = -Icges(6,3) * t221 + t243 * t220;
t201 = t231 * (rSges(4,1) * t223 - rSges(4,2) * t222) + t251;
t200 = -t248 + (-rSges(4,1) * t222 - rSges(4,2) * t223 - t255) * t231;
t199 = t226 * (rSges(5,1) * t221 - rSges(5,2) * t220) + t247;
t198 = -t226 * (rSges(5,1) * t220 + rSges(5,2) * t221) + t239;
t197 = qJD(3) + (t208 * t220 + t209 * t221) * qJD(5);
t196 = -t217 * t250 + (pkin(4) * t221 + pkin(8) * t220 + t209) * t226 + t247;
t195 = -t217 * t249 + (-pkin(4) * t220 + pkin(8) * t221 - t208) * t226 + t239;
t1 = m(3) * (t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(4) * (t200 ^ 2 + t201 ^ 2 + t237) / 0.2e1 + m(5) * (t198 ^ 2 + t199 ^ 2 + t237) / 0.2e1 + t226 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + ((t220 * t213 + t240 * t221) * t226 + (t220 ^ 2 * t203 + (t242 * t221 + (-t202 + t241) * t220) * t221) * qJD(5)) * t250 / 0.2e1 - ((-t221 * t213 + t240 * t220) * t226 + (t221 ^ 2 * t202 + (t241 * t220 + (-t203 + t242) * t221) * t220) * qJD(5)) * t249 / 0.2e1 + t226 * ((t235 * t214 + t233 * t215) * t226 + ((t205 * t235 + t207 * t233) * t220 - (t204 * t235 + t206 * t233) * t221) * qJD(5)) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t231 ^ 2 / 0.2e1 + (m(2) * (t218 ^ 2 + t219 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
