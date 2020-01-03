% Calculate kinetic energy for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:35
% EndTime: 2019-12-31 16:41:35
% DurationCPUTime: 0.31s
% Computational Cost: add. (211->85), mult. (293->139), div. (0->0), fcn. (226->6), ass. (0->48)
t207 = sin(qJ(1));
t208 = cos(qJ(1));
t225 = qJD(1) * t208 * qJ(3) + qJD(3) * t207;
t204 = sin(pkin(6));
t223 = pkin(3) * t204;
t203 = pkin(6) + qJ(4);
t198 = sin(t203);
t222 = Icges(5,4) * t198;
t199 = cos(t203);
t221 = Icges(5,4) * t199;
t202 = qJD(2) * t207;
t220 = qJD(3) * t208 + t202;
t219 = qJD(4) * t207;
t193 = qJD(1) * (pkin(1) * t208 + qJ(2) * t207);
t218 = -qJD(2) * t208 + t193;
t205 = cos(pkin(6));
t217 = rSges(4,1) * t204 + rSges(4,2) * t205;
t216 = rSges(5,1) * t198 + rSges(5,2) * t199;
t215 = Icges(5,1) * t198 + t221;
t214 = Icges(5,2) * t199 + t222;
t213 = Icges(5,5) * t198 + Icges(5,6) * t199;
t183 = Icges(5,6) * t208 + t214 * t207;
t185 = Icges(5,5) * t208 + t215 * t207;
t212 = -t183 * t199 - t185 * t198;
t184 = Icges(5,6) * t207 - t214 * t208;
t186 = Icges(5,5) * t207 - t215 * t208;
t211 = t184 * t199 + t186 * t198;
t190 = -Icges(5,2) * t198 + t221;
t191 = Icges(5,1) * t199 - t222;
t210 = t190 * t199 + t191 * t198;
t206 = -pkin(5) - qJ(3);
t196 = rSges(2,1) * t208 - rSges(2,2) * t207;
t195 = rSges(2,1) * t207 + rSges(2,2) * t208;
t194 = pkin(1) * t207 - qJ(2) * t208;
t192 = rSges(5,1) * t199 - rSges(5,2) * t198;
t189 = Icges(5,5) * t199 - Icges(5,6) * t198;
t188 = rSges(5,3) * t207 - t216 * t208;
t187 = rSges(5,3) * t208 + t216 * t207;
t182 = Icges(5,3) * t207 - t213 * t208;
t181 = Icges(5,3) * t208 + t213 * t207;
t180 = qJD(1) * (-rSges(3,2) * t208 + rSges(3,3) * t207) + t218;
t179 = t202 + (rSges(3,2) * t207 + rSges(3,3) * t208 - t194) * qJD(1);
t178 = qJD(1) * (rSges(4,3) * t208 + t217 * t207) + t218 + t225;
t177 = (-t194 + t217 * t208 + (-rSges(4,3) - qJ(3)) * t207) * qJD(1) + t220;
t176 = (-t187 * t207 + t188 * t208) * qJD(4);
t175 = t193 + (t207 * t223 + t187) * qJD(1) + (-qJD(2) + qJD(1) * (-qJ(3) - t206) - qJD(4) * t192) * t208 + t225;
t174 = t192 * t219 + (t206 * t207 + t208 * t223 - t188 - t194) * qJD(1) + t220;
t1 = m(3) * (t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(4) * (t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(5) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + qJD(4) * t208 * ((t208 * t189 + t210 * t207) * qJD(1) + (t208 ^ 2 * t181 + (t211 * t207 + (t182 - t212) * t208) * t207) * qJD(4)) / 0.2e1 + ((t207 * t189 - t210 * t208) * qJD(1) + (t207 ^ 2 * t182 + (t212 * t208 + (t181 - t211) * t207) * t208) * qJD(4)) * t219 / 0.2e1 + qJD(1) * ((-t198 * t190 + t199 * t191) * qJD(1) + ((-t183 * t198 + t185 * t199) * t208 + (-t184 * t198 + t186 * t199) * t207) * qJD(4)) / 0.2e1 + (m(2) * (t195 ^ 2 + t196 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t205 ^ 2 + (-0.2e1 * Icges(4,4) * t205 + Icges(4,2) * t204) * t204) * qJD(1) ^ 2 / 0.2e1;
T = t1;
