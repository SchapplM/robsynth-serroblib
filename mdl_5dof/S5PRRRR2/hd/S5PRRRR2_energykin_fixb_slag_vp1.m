% Calculate kinetic energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:10
% EndTime: 2019-07-18 13:30:10
% DurationCPUTime: 0.37s
% Computational Cost: add. (340->57), mult. (265->106), div. (0->0), fcn. (192->8), ass. (0->42)
t193 = qJ(2) + qJ(3);
t191 = qJ(4) + t193;
t185 = sin(t191);
t219 = t185 ^ 2;
t186 = cos(t191);
t218 = t186 ^ 2;
t194 = sin(qJ(5));
t196 = cos(qJ(5));
t217 = Icges(6,5) * t194 + Icges(6,6) * t196;
t192 = qJD(2) + qJD(3);
t188 = qJD(4) + t192;
t216 = t217 * t188;
t215 = t185 * t186;
t214 = pkin(3) * t192;
t213 = pkin(2) * qJD(2);
t197 = cos(qJ(2));
t187 = t197 * t213;
t190 = cos(t193);
t210 = t190 * t214 + t187;
t195 = sin(qJ(2));
t209 = t195 * t213;
t208 = rSges(6,1) * t196 - rSges(6,2) * t194;
t207 = pkin(6) * t188 - qJD(5) * (t194 * rSges(6,1) + t196 * rSges(6,2));
t204 = Icges(6,5) * t196 - Icges(6,6) * t194;
t189 = sin(t193);
t200 = -t189 * t214 - t209;
t199 = qJD(1) ^ 2;
t198 = qJD(2) ^ 2;
t184 = t197 * rSges(3,1) - t195 * rSges(3,2);
t183 = t195 * rSges(3,1) + t197 * rSges(3,2);
t177 = t187 + t192 * (t190 * rSges(4,1) - t189 * rSges(4,2));
t176 = -t209 - t192 * (t189 * rSges(4,1) + t190 * rSges(4,2));
t175 = t185 * rSges(6,3) + t208 * t186;
t174 = -t186 * rSges(6,3) + t208 * t185;
t169 = Icges(6,3) * t185 + t204 * t186;
t168 = -Icges(6,3) * t186 + t204 * t185;
t167 = t188 * (t186 * rSges(5,1) - t185 * rSges(5,2)) + t210;
t166 = -t188 * (t185 * rSges(5,1) + t186 * rSges(5,2)) + t200;
t165 = qJD(1) + (t174 * t185 + t175 * t186) * qJD(5);
t164 = t188 * t175 + t207 * t185 + t210;
t163 = -t188 * t174 + t207 * t186 + t200;
t1 = m(2) * t199 / 0.2e1 + m(3) * (t199 + (t183 ^ 2 + t184 ^ 2) * t198) / 0.2e1 + t198 * Icges(3,3) / 0.2e1 + m(4) * (t176 ^ 2 + t177 ^ 2 + t199) / 0.2e1 + t192 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t166 ^ 2 + t167 ^ 2 + t199) / 0.2e1 + t188 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + qJD(5) * t185 * (t185 * t216 + (-t168 * t215 + t219 * t169) * qJD(5)) / 0.2e1 - qJD(5) * t186 * (-t186 * t216 + (t218 * t168 - t169 * t215) * qJD(5)) / 0.2e1 + t188 * ((Icges(6,2) * t196 ^ 2 + (Icges(6,1) * t194 + 0.2e1 * Icges(6,4) * t196) * t194) * t188 + (t218 + t219) * t217 * qJD(5)) / 0.2e1;
T  = t1;
