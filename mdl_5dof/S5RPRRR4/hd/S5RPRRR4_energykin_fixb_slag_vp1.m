% Calculate kinetic energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:27
% EndTime: 2022-01-23 09:34:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (474->87), mult. (293->144), div. (0->0), fcn. (206->10), ass. (0->58)
t232 = sin(qJ(1));
t253 = pkin(1) * t232;
t229 = qJD(1) + qJD(3);
t252 = pkin(3) * t229;
t231 = sin(qJ(5));
t251 = Icges(6,4) * t231;
t233 = cos(qJ(5));
t250 = Icges(6,4) * t233;
t234 = cos(qJ(1));
t224 = qJD(1) * t234 * pkin(1);
t230 = qJ(1) + pkin(9);
t226 = cos(t230);
t249 = qJD(1) * pkin(2) * t226 + t224;
t228 = qJ(3) + t230;
t223 = qJ(4) + t228;
t219 = sin(t223);
t248 = qJD(5) * t219;
t220 = cos(t223);
t247 = qJD(5) * t220;
t222 = cos(t228);
t246 = t222 * t252 + t249;
t245 = rSges(6,1) * t233 - rSges(6,2) * t231;
t244 = Icges(6,1) * t233 - t251;
t243 = -Icges(6,2) * t231 + t250;
t242 = Icges(6,5) * t233 - Icges(6,6) * t231;
t203 = -Icges(6,6) * t220 + t243 * t219;
t205 = -Icges(6,5) * t220 + t244 * t219;
t241 = t203 * t231 - t205 * t233;
t204 = Icges(6,6) * t219 + t243 * t220;
t206 = Icges(6,5) * t219 + t244 * t220;
t240 = -t204 * t231 + t206 * t233;
t213 = Icges(6,2) * t233 + t251;
t214 = Icges(6,1) * t231 + t250;
t239 = -t213 * t231 + t214 * t233;
t225 = sin(t230);
t238 = (-pkin(2) * t225 - t253) * qJD(1);
t221 = sin(t228);
t237 = -t221 * t252 + t238;
t235 = qJD(2) ^ 2;
t227 = qJD(4) + t229;
t217 = rSges(2,1) * t234 - rSges(2,2) * t232;
t216 = rSges(2,1) * t232 + rSges(2,2) * t234;
t215 = rSges(6,1) * t231 + rSges(6,2) * t233;
t212 = Icges(6,5) * t231 + Icges(6,6) * t233;
t210 = t224 + qJD(1) * (rSges(3,1) * t226 - rSges(3,2) * t225);
t209 = (-rSges(3,1) * t225 - rSges(3,2) * t226 - t253) * qJD(1);
t208 = rSges(6,3) * t219 + t245 * t220;
t207 = -rSges(6,3) * t220 + t245 * t219;
t202 = Icges(6,3) * t219 + t242 * t220;
t201 = -Icges(6,3) * t220 + t242 * t219;
t200 = t229 * (rSges(4,1) * t222 - rSges(4,2) * t221) + t249;
t199 = -t229 * (rSges(4,1) * t221 + rSges(4,2) * t222) + t238;
t198 = t227 * (rSges(5,1) * t220 - rSges(5,2) * t219) + t246;
t197 = -t227 * (rSges(5,1) * t219 + rSges(5,2) * t220) + t237;
t196 = qJD(2) + (t207 * t219 + t208 * t220) * qJD(5);
t195 = -t215 * t248 + (pkin(4) * t220 + pkin(8) * t219 + t208) * t227 + t246;
t194 = -t215 * t247 + (-pkin(4) * t219 + pkin(8) * t220 - t207) * t227 + t237;
t1 = m(3) * (t209 ^ 2 + t210 ^ 2 + t235) / 0.2e1 + m(4) * (t199 ^ 2 + t200 ^ 2 + t235) / 0.2e1 + t229 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t197 ^ 2 + t198 ^ 2 + t235) / 0.2e1 + t227 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + ((t219 * t212 + t239 * t220) * t227 + (t219 ^ 2 * t202 + (t241 * t220 + (-t201 + t240) * t219) * t220) * qJD(5)) * t248 / 0.2e1 - ((-t220 * t212 + t239 * t219) * t227 + (t220 ^ 2 * t201 + (t240 * t219 + (-t202 + t241) * t220) * t219) * qJD(5)) * t247 / 0.2e1 + t227 * ((t213 * t233 + t214 * t231) * t227 + ((t204 * t233 + t206 * t231) * t219 - (t203 * t233 + t205 * t231) * t220) * qJD(5)) / 0.2e1 + (m(2) * (t216 ^ 2 + t217 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
