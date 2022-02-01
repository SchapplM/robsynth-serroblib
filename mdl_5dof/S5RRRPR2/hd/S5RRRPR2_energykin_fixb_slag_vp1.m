% Calculate kinetic energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:31
% EndTime: 2022-01-20 11:30:31
% DurationCPUTime: 0.28s
% Computational Cost: add. (488->86), mult. (291->144), div. (0->0), fcn. (206->10), ass. (0->58)
t230 = qJD(1) + qJD(2);
t254 = pkin(2) * t230;
t231 = qJ(1) + qJ(2);
t229 = qJ(3) + t231;
t222 = sin(t229);
t253 = pkin(3) * t222;
t252 = pkin(1) * qJD(1);
t232 = sin(qJ(5));
t251 = Icges(6,4) * t232;
t234 = cos(qJ(5));
t250 = Icges(6,4) * t234;
t235 = cos(qJ(1));
t224 = t235 * t252;
t228 = cos(t231);
t249 = t228 * t254 + t224;
t221 = pkin(9) + t229;
t219 = sin(t221);
t248 = qJD(5) * t219;
t220 = cos(t221);
t247 = qJD(5) * t220;
t233 = sin(qJ(1));
t246 = t233 * t252;
t223 = cos(t229);
t226 = qJD(3) + t230;
t245 = t226 * pkin(3) * t223 + t249;
t244 = rSges(6,1) * t234 - rSges(6,2) * t232;
t243 = Icges(6,1) * t234 - t251;
t242 = -Icges(6,2) * t232 + t250;
t241 = Icges(6,5) * t234 - Icges(6,6) * t232;
t203 = -Icges(6,6) * t220 + t242 * t219;
t205 = -Icges(6,5) * t220 + t243 * t219;
t240 = t203 * t232 - t205 * t234;
t204 = Icges(6,6) * t219 + t242 * t220;
t206 = Icges(6,5) * t219 + t243 * t220;
t239 = -t204 * t232 + t206 * t234;
t213 = Icges(6,2) * t234 + t251;
t214 = Icges(6,1) * t232 + t250;
t238 = -t213 * t232 + t214 * t234;
t227 = sin(t231);
t237 = -t227 * t254 - t246;
t218 = t235 * rSges(2,1) - t233 * rSges(2,2);
t217 = t233 * rSges(2,1) + t235 * rSges(2,2);
t216 = t232 * rSges(6,1) + t234 * rSges(6,2);
t212 = Icges(6,5) * t232 + Icges(6,6) * t234;
t210 = t224 + t230 * (t228 * rSges(3,1) - t227 * rSges(3,2));
t209 = -t246 - t230 * (t227 * rSges(3,1) + t228 * rSges(3,2));
t208 = t219 * rSges(6,3) + t244 * t220;
t207 = -t220 * rSges(6,3) + t244 * t219;
t202 = Icges(6,3) * t219 + t241 * t220;
t201 = -Icges(6,3) * t220 + t241 * t219;
t200 = t226 * (t223 * rSges(4,1) - t222 * rSges(4,2)) + t249;
t199 = -t226 * (t222 * rSges(4,1) + t223 * rSges(4,2)) + t237;
t198 = t226 * (t220 * rSges(5,1) - t219 * rSges(5,2)) + t245;
t197 = (-t219 * rSges(5,1) - t220 * rSges(5,2) - t253) * t226 + t237;
t196 = qJD(4) + (t207 * t219 + t208 * t220) * qJD(5);
t195 = -t216 * t248 + (t220 * pkin(4) + t219 * pkin(8) + t208) * t226 + t245;
t194 = -t216 * t247 + (-t219 * pkin(4) + t220 * pkin(8) - t207 - t253) * t226 + t237;
t1 = m(3) * (t209 ^ 2 + t210 ^ 2) / 0.2e1 + t230 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t199 ^ 2 + t200 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(6) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + ((t219 * t212 + t238 * t220) * t226 + (t219 ^ 2 * t202 + (t240 * t220 + (-t201 + t239) * t219) * t220) * qJD(5)) * t248 / 0.2e1 - ((-t220 * t212 + t238 * t219) * t226 + (t220 ^ 2 * t201 + (t239 * t219 + (-t202 + t240) * t220) * t219) * qJD(5)) * t247 / 0.2e1 + t226 * ((t234 * t213 + t232 * t214) * t226 + ((t234 * t204 + t232 * t206) * t219 - (t234 * t203 + t232 * t205) * t220) * qJD(5)) / 0.2e1 + (Icges(4,3) + Icges(5,3)) * t226 ^ 2 / 0.2e1 + (m(2) * (t217 ^ 2 + t218 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
