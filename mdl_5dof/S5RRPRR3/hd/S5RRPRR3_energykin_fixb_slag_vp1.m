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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:00:20
% EndTime: 2020-01-03 12:00:21
% DurationCPUTime: 0.26s
% Computational Cost: add. (482->86), mult. (292->144), div. (0->0), fcn. (206->10), ass. (0->60)
t231 = qJD(1) + qJD(2);
t256 = pkin(2) * t231;
t255 = pkin(3) * t231;
t254 = pkin(1) * qJD(1);
t233 = sin(qJ(5));
t253 = Icges(6,4) * t233;
t235 = cos(qJ(5));
t252 = Icges(6,4) * t235;
t232 = qJ(1) + qJ(2);
t234 = sin(qJ(1));
t224 = t234 * t254;
t228 = sin(t232);
t251 = t228 * t256 + t224;
t236 = cos(qJ(1));
t225 = t236 * t254;
t229 = cos(t232);
t250 = t229 * t256 + t225;
t227 = pkin(9) + t232;
t223 = qJ(4) + t227;
t219 = sin(t223);
t249 = qJD(5) * t219;
t220 = cos(t223);
t248 = qJD(5) * t220;
t221 = sin(t227);
t247 = t221 * t255 + t251;
t222 = cos(t227);
t246 = t222 * t255 + t250;
t245 = rSges(6,1) * t235 - rSges(6,2) * t233;
t244 = Icges(6,1) * t235 - t253;
t243 = -Icges(6,2) * t233 + t252;
t242 = Icges(6,5) * t235 - Icges(6,6) * t233;
t201 = -Icges(6,6) * t220 + t243 * t219;
t203 = -Icges(6,5) * t220 + t244 * t219;
t241 = -t201 * t233 + t203 * t235;
t202 = -Icges(6,6) * t219 - t243 * t220;
t204 = -Icges(6,5) * t219 - t244 * t220;
t240 = t202 * t233 - t204 * t235;
t212 = Icges(6,2) * t235 + t253;
t213 = Icges(6,1) * t233 + t252;
t239 = t212 * t233 - t213 * t235;
t237 = qJD(3) ^ 2;
t226 = qJD(4) + t231;
t218 = -t236 * rSges(2,1) + t234 * rSges(2,2);
t217 = t234 * rSges(2,1) + t236 * rSges(2,2);
t216 = t233 * rSges(6,1) + t235 * rSges(6,2);
t211 = Icges(6,5) * t233 + Icges(6,6) * t235;
t208 = t225 - t231 * (-t229 * rSges(3,1) + t228 * rSges(3,2));
t207 = t224 + t231 * (t228 * rSges(3,1) + t229 * rSges(3,2));
t206 = -t219 * rSges(6,3) - t245 * t220;
t205 = -t220 * rSges(6,3) + t245 * t219;
t200 = -Icges(6,3) * t219 - t242 * t220;
t199 = -Icges(6,3) * t220 + t242 * t219;
t198 = -t231 * (-t222 * rSges(4,1) + t221 * rSges(4,2)) + t250;
t197 = t231 * (t221 * rSges(4,1) + t222 * rSges(4,2)) + t251;
t196 = -t226 * (-t220 * rSges(5,1) + t219 * rSges(5,2)) + t246;
t195 = t226 * (t219 * rSges(5,1) + t220 * rSges(5,2)) + t247;
t194 = qJD(3) + (t205 * t219 - t206 * t220) * qJD(5);
t193 = -t216 * t249 + (t220 * pkin(4) + t219 * pkin(8) - t206) * t226 + t246;
t192 = t216 * t248 + (t219 * pkin(4) - t220 * pkin(8) + t205) * t226 + t247;
t1 = m(3) * (t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(4) * (t197 ^ 2 + t198 ^ 2 + t237) / 0.2e1 + m(5) * (t195 ^ 2 + t196 ^ 2 + t237) / 0.2e1 + t226 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + t226 * ((t235 * t212 + t233 * t213) * t226 + (-(t235 * t201 + t233 * t203) * t220 - (t235 * t202 + t233 * t204) * t219) * qJD(5)) / 0.2e1 - ((-t220 * t211 - t239 * t219) * t226 + (t220 ^ 2 * t199 + (t240 * t219 + (t200 - t241) * t220) * t219) * qJD(5)) * t248 / 0.2e1 - ((-t219 * t211 + t239 * t220) * t226 + (t219 ^ 2 * t200 + (t241 * t220 + (t199 - t240) * t219) * t220) * qJD(5)) * t249 / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t231 ^ 2 / 0.2e1 + (m(2) * (t217 ^ 2 + t218 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
