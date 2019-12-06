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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:17
% EndTime: 2019-12-05 18:30:17
% DurationCPUTime: 0.28s
% Computational Cost: add. (482->88), mult. (292->144), div. (0->0), fcn. (206->10), ass. (0->58)
t224 = qJ(1) + qJ(2);
t220 = sin(t224);
t248 = pkin(2) * t220;
t221 = cos(t224);
t247 = pkin(2) * t221;
t246 = pkin(1) * qJD(1);
t225 = sin(qJ(5));
t245 = Icges(6,4) * t225;
t227 = cos(qJ(5));
t244 = Icges(6,4) * t227;
t219 = pkin(9) + t224;
t217 = qJ(4) + t219;
t213 = sin(t217);
t243 = qJD(5) * t213;
t214 = cos(t217);
t242 = qJD(5) * t214;
t223 = qJD(1) + qJD(2);
t226 = sin(qJ(1));
t241 = t226 * t246;
t228 = cos(qJ(1));
t240 = t228 * t246;
t239 = rSges(6,1) * t227 - rSges(6,2) * t225;
t238 = Icges(6,1) * t227 - t245;
t237 = -Icges(6,2) * t225 + t244;
t236 = Icges(6,5) * t227 - Icges(6,6) * t225;
t199 = Icges(6,6) * t214 - t237 * t213;
t201 = Icges(6,5) * t214 - t238 * t213;
t235 = -t199 * t225 + t201 * t227;
t200 = Icges(6,6) * t213 + t237 * t214;
t202 = Icges(6,5) * t213 + t238 * t214;
t234 = t200 * t225 - t202 * t227;
t208 = Icges(6,2) * t227 + t245;
t209 = Icges(6,1) * t225 + t244;
t233 = t208 * t225 - t209 * t227;
t215 = sin(t219);
t232 = (-pkin(3) * t215 - t248) * t223 - t241;
t216 = cos(t219);
t231 = (-pkin(3) * t216 - t247) * t223 - t240;
t229 = qJD(3) ^ 2;
t218 = qJD(4) + t223;
t212 = t228 * rSges(2,1) - t226 * rSges(2,2);
t211 = -t226 * rSges(2,1) - t228 * rSges(2,2);
t210 = t225 * rSges(6,1) + t227 * rSges(6,2);
t207 = Icges(6,5) * t225 + Icges(6,6) * t227;
t206 = -t240 - t223 * (t221 * rSges(3,1) - t220 * rSges(3,2));
t205 = -t241 + t223 * (-t220 * rSges(3,1) - t221 * rSges(3,2));
t204 = t213 * rSges(6,3) + t239 * t214;
t203 = t214 * rSges(6,3) - t239 * t213;
t198 = Icges(6,3) * t213 + t236 * t214;
t197 = Icges(6,3) * t214 - t236 * t213;
t196 = -t240 + (-t216 * rSges(4,1) + t215 * rSges(4,2) - t247) * t223;
t195 = -t241 + (-t215 * rSges(4,1) - t216 * rSges(4,2) - t248) * t223;
t194 = -t218 * (t214 * rSges(5,1) - t213 * rSges(5,2)) + t231;
t193 = t218 * (-t213 * rSges(5,1) - t214 * rSges(5,2)) + t232;
t192 = qJD(3) + (-t203 * t213 + t204 * t214) * qJD(5);
t191 = t210 * t243 + (-t214 * pkin(4) - t213 * pkin(8) - t204) * t218 + t231;
t190 = -t210 * t242 + (-t213 * pkin(4) + t214 * pkin(8) + t203) * t218 + t232;
t1 = m(3) * (t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(4) * (t195 ^ 2 + t196 ^ 2 + t229) / 0.2e1 + m(5) * (t193 ^ 2 + t194 ^ 2 + t229) / 0.2e1 + t218 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + t218 * ((t227 * t208 + t225 * t209) * t218 + ((t227 * t199 + t225 * t201) * t214 + (t227 * t200 + t225 * t202) * t213) * qJD(5)) / 0.2e1 + ((t214 * t207 + t233 * t213) * t218 + (t214 ^ 2 * t197 + (t234 * t213 + (t198 - t235) * t214) * t213) * qJD(5)) * t242 / 0.2e1 + ((t213 * t207 - t233 * t214) * t218 + (t213 ^ 2 * t198 + (t235 * t214 + (t197 - t234) * t213) * t214) * qJD(5)) * t243 / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t223 ^ 2 / 0.2e1 + (m(2) * (t211 ^ 2 + t212 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
