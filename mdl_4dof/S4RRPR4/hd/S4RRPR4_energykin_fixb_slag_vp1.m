% Calculate kinetic energy for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:25
% EndTime: 2019-12-31 17:02:26
% DurationCPUTime: 0.30s
% Computational Cost: add. (345->81), mult. (287->142), div. (0->0), fcn. (220->8), ass. (0->51)
t198 = cos(pkin(7));
t219 = t198 * pkin(3);
t218 = pkin(1) * qJD(1);
t194 = pkin(7) + qJ(4);
t189 = sin(t194);
t217 = Icges(5,4) * t189;
t190 = cos(t194);
t216 = Icges(5,4) * t190;
t201 = cos(qJ(1));
t188 = t201 * t218;
t196 = qJ(1) + qJ(2);
t191 = sin(t196);
t192 = cos(t196);
t195 = qJD(1) + qJD(2);
t214 = t195 * (t192 * pkin(2) + t191 * qJ(3)) + t188;
t213 = qJD(4) * t192;
t200 = sin(qJ(1));
t212 = t200 * t218;
t211 = qJD(3) * t191 - t212;
t197 = sin(pkin(7));
t210 = rSges(4,1) * t198 - rSges(4,2) * t197;
t209 = rSges(5,1) * t190 - rSges(5,2) * t189;
t208 = Icges(5,1) * t190 - t217;
t207 = -Icges(5,2) * t189 + t216;
t206 = Icges(5,5) * t190 - Icges(5,6) * t189;
t170 = -Icges(5,6) * t192 + t207 * t191;
t172 = -Icges(5,5) * t192 + t208 * t191;
t205 = t170 * t189 - t172 * t190;
t171 = Icges(5,6) * t191 + t192 * t207;
t173 = Icges(5,5) * t191 + t192 * t208;
t204 = -t171 * t189 + t173 * t190;
t180 = Icges(5,2) * t190 + t217;
t181 = Icges(5,1) * t189 + t216;
t203 = -t180 * t189 + t181 * t190;
t185 = t201 * rSges(2,1) - t200 * rSges(2,2);
t184 = t200 * rSges(2,1) + t201 * rSges(2,2);
t183 = t191 * pkin(2) - t192 * qJ(3);
t182 = t189 * rSges(5,1) + t190 * rSges(5,2);
t179 = Icges(5,5) * t189 + Icges(5,6) * t190;
t177 = t188 + t195 * (t192 * rSges(3,1) - t191 * rSges(3,2));
t176 = -t212 - t195 * (t191 * rSges(3,1) + t192 * rSges(3,2));
t175 = t191 * rSges(5,3) + t192 * t209;
t174 = -t192 * rSges(5,3) + t209 * t191;
t169 = Icges(5,3) * t191 + t192 * t206;
t168 = -Icges(5,3) * t192 + t206 * t191;
t167 = t195 * t191 * rSges(4,3) + (t195 * t210 - qJD(3)) * t192 + t214;
t166 = (t192 * rSges(4,3) - t210 * t191 - t183) * t195 + t211;
t165 = (t174 * t191 + t175 * t192) * qJD(4);
t164 = t195 * t175 + (t195 * t219 - qJD(3)) * t192 + (pkin(6) * t195 - qJD(4) * t182) * t191 + t214;
t163 = -t182 * t213 + (pkin(6) * t192 - t191 * t219 - t174 - t183) * t195 + t211;
t1 = m(3) * (t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(4) * (t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(5) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + qJD(4) * t191 * ((t191 * t179 + t192 * t203) * t195 + (t191 ^ 2 * t169 + (t205 * t192 + (-t168 + t204) * t191) * t192) * qJD(4)) / 0.2e1 - ((-t192 * t179 + t203 * t191) * t195 + (t192 ^ 2 * t168 + (t204 * t191 + (-t169 + t205) * t192) * t191) * qJD(4)) * t213 / 0.2e1 + t195 * ((t190 * t180 + t189 * t181) * t195 + ((t190 * t171 + t189 * t173) * t191 - (t190 * t170 + t189 * t172) * t192) * qJD(4)) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t198 ^ 2 + (Icges(4,1) * t197 + 0.2e1 * Icges(4,4) * t198) * t197) * t195 ^ 2 / 0.2e1 + (m(2) * (t184 ^ 2 + t185 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
