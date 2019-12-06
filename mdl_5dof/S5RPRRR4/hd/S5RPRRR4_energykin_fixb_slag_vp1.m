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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:14:18
% EndTime: 2019-12-05 18:14:18
% DurationCPUTime: 0.29s
% Computational Cost: add. (474->87), mult. (293->144), div. (0->0), fcn. (206->10), ass. (0->58)
t221 = qJD(1) + qJD(3);
t246 = pkin(3) * t221;
t224 = sin(qJ(1));
t245 = t224 * pkin(1);
t226 = cos(qJ(1));
t244 = t226 * pkin(1);
t223 = sin(qJ(5));
t243 = Icges(6,4) * t223;
t225 = cos(qJ(5));
t242 = Icges(6,4) * t225;
t222 = qJ(1) + pkin(9);
t220 = qJ(3) + t222;
t216 = qJ(4) + t220;
t212 = sin(t216);
t241 = qJD(5) * t212;
t213 = cos(t216);
t240 = qJD(5) * t213;
t239 = rSges(6,1) * t225 - rSges(6,2) * t223;
t238 = Icges(6,1) * t225 - t243;
t237 = -Icges(6,2) * t223 + t242;
t236 = Icges(6,5) * t225 - Icges(6,6) * t223;
t198 = Icges(6,6) * t213 - t237 * t212;
t200 = Icges(6,5) * t213 - t238 * t212;
t235 = -t198 * t223 + t200 * t225;
t199 = Icges(6,6) * t212 + t237 * t213;
t201 = Icges(6,5) * t212 + t238 * t213;
t234 = t199 * t223 - t201 * t225;
t207 = Icges(6,2) * t225 + t243;
t208 = Icges(6,1) * t223 + t242;
t233 = t207 * t223 - t208 * t225;
t217 = sin(t222);
t232 = (-pkin(2) * t217 - t245) * qJD(1);
t218 = cos(t222);
t231 = (-pkin(2) * t218 - t244) * qJD(1);
t214 = sin(t220);
t230 = -t214 * t246 + t232;
t215 = cos(t220);
t229 = -t215 * t246 + t231;
t227 = qJD(2) ^ 2;
t219 = qJD(4) + t221;
t211 = t226 * rSges(2,1) - t224 * rSges(2,2);
t210 = -t224 * rSges(2,1) - t226 * rSges(2,2);
t209 = t223 * rSges(6,1) + t225 * rSges(6,2);
t206 = Icges(6,5) * t223 + Icges(6,6) * t225;
t205 = (-t218 * rSges(3,1) + t217 * rSges(3,2) - t244) * qJD(1);
t204 = (-t217 * rSges(3,1) - t218 * rSges(3,2) - t245) * qJD(1);
t203 = t212 * rSges(6,3) + t239 * t213;
t202 = t213 * rSges(6,3) - t239 * t212;
t197 = Icges(6,3) * t212 + t236 * t213;
t196 = Icges(6,3) * t213 - t236 * t212;
t195 = -t221 * (t215 * rSges(4,1) - t214 * rSges(4,2)) + t231;
t194 = t221 * (-t214 * rSges(4,1) - t215 * rSges(4,2)) + t232;
t193 = -t219 * (t213 * rSges(5,1) - t212 * rSges(5,2)) + t229;
t192 = t219 * (-t212 * rSges(5,1) - t213 * rSges(5,2)) + t230;
t191 = qJD(2) + (-t202 * t212 + t203 * t213) * qJD(5);
t190 = t209 * t241 + (-t213 * pkin(4) - t212 * pkin(8) - t203) * t219 + t229;
t189 = -t209 * t240 + (-t212 * pkin(4) + t213 * pkin(8) + t202) * t219 + t230;
t1 = m(3) * (t204 ^ 2 + t205 ^ 2 + t227) / 0.2e1 + m(4) * (t194 ^ 2 + t195 ^ 2 + t227) / 0.2e1 + t221 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t192 ^ 2 + t193 ^ 2 + t227) / 0.2e1 + t219 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + t219 * ((t225 * t207 + t223 * t208) * t219 + ((t225 * t198 + t223 * t200) * t213 + (t225 * t199 + t223 * t201) * t212) * qJD(5)) / 0.2e1 + ((t213 * t206 + t233 * t212) * t219 + (t213 ^ 2 * t196 + (t234 * t212 + (t197 - t235) * t213) * t212) * qJD(5)) * t240 / 0.2e1 + ((t212 * t206 - t233 * t213) * t219 + (t212 ^ 2 * t197 + (t235 * t213 + (t196 - t234) * t212) * t213) * qJD(5)) * t241 / 0.2e1 + (m(2) * (t210 ^ 2 + t211 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
