% Calculate kinetic energy for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:03
% EndTime: 2019-12-05 17:06:03
% DurationCPUTime: 0.23s
% Computational Cost: add. (462->79), mult. (267->135), div. (0->0), fcn. (194->8), ass. (0->54)
t227 = qJD(2) + qJD(3);
t247 = pkin(3) * t227;
t246 = pkin(2) * qJD(2);
t228 = sin(qJ(5));
t245 = Icges(6,4) * t228;
t229 = cos(qJ(5));
t244 = Icges(6,4) * t229;
t226 = pkin(9) + qJ(2);
t223 = cos(t226);
t216 = t223 * t246;
t225 = qJ(3) + t226;
t220 = cos(t225);
t243 = t220 * t247 + t216;
t221 = qJ(4) + t225;
t217 = sin(t221);
t242 = qJD(5) * t217;
t218 = cos(t221);
t241 = qJD(5) * t218;
t222 = sin(t226);
t240 = t222 * t246;
t239 = rSges(6,1) * t229 - rSges(6,2) * t228;
t238 = Icges(6,1) * t229 - t245;
t237 = -Icges(6,2) * t228 + t244;
t236 = Icges(6,5) * t229 - Icges(6,6) * t228;
t201 = -Icges(6,6) * t218 + t237 * t217;
t203 = -Icges(6,5) * t218 + t238 * t217;
t235 = t201 * t228 - t203 * t229;
t202 = Icges(6,6) * t217 + t237 * t218;
t204 = Icges(6,5) * t217 + t238 * t218;
t234 = -t202 * t228 + t204 * t229;
t213 = Icges(6,2) * t229 + t245;
t214 = Icges(6,1) * t228 + t244;
t233 = -t213 * t228 + t214 * t229;
t219 = sin(t225);
t232 = -t219 * t247 - t240;
t231 = qJD(1) ^ 2;
t230 = qJD(2) ^ 2;
t224 = qJD(4) + t227;
t215 = rSges(6,1) * t228 + rSges(6,2) * t229;
t212 = Icges(6,5) * t228 + Icges(6,6) * t229;
t210 = rSges(3,1) * t223 - rSges(3,2) * t222;
t209 = rSges(3,1) * t222 + rSges(3,2) * t223;
t208 = t216 + t227 * (rSges(4,1) * t220 - rSges(4,2) * t219);
t207 = -t240 - t227 * (rSges(4,1) * t219 + rSges(4,2) * t220);
t206 = rSges(6,3) * t217 + t239 * t218;
t205 = -rSges(6,3) * t218 + t239 * t217;
t200 = Icges(6,3) * t217 + t236 * t218;
t199 = -Icges(6,3) * t218 + t236 * t217;
t198 = t224 * (rSges(5,1) * t218 - rSges(5,2) * t217) + t243;
t197 = -t224 * (rSges(5,1) * t217 + rSges(5,2) * t218) + t232;
t196 = qJD(1) + (t205 * t217 + t206 * t218) * qJD(5);
t195 = -t215 * t242 + (pkin(4) * t218 + pkin(8) * t217 + t206) * t224 + t243;
t194 = -t215 * t241 + (-pkin(4) * t217 + pkin(8) * t218 - t205) * t224 + t232;
t1 = m(2) * t231 / 0.2e1 + m(3) * (t231 + (t209 ^ 2 + t210 ^ 2) * t230) / 0.2e1 + t230 * Icges(3,3) / 0.2e1 + m(4) * (t207 ^ 2 + t208 ^ 2 + t231) / 0.2e1 + t227 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t197 ^ 2 + t198 ^ 2 + t231) / 0.2e1 + t224 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + ((t217 * t212 + t233 * t218) * t224 + (t217 ^ 2 * t200 + (t235 * t218 + (-t199 + t234) * t217) * t218) * qJD(5)) * t242 / 0.2e1 - ((-t218 * t212 + t233 * t217) * t224 + (t218 ^ 2 * t199 + (t234 * t217 + (-t200 + t235) * t218) * t217) * qJD(5)) * t241 / 0.2e1 + t224 * ((t229 * t213 + t228 * t214) * t224 + ((t202 * t229 + t204 * t228) * t217 - (t201 * t229 + t203 * t228) * t218) * qJD(5)) / 0.2e1;
T = t1;
