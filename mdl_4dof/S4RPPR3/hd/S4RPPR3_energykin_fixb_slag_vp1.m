% Calculate kinetic energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:49
% EndTime: 2019-12-31 16:37:49
% DurationCPUTime: 0.26s
% Computational Cost: add. (330->83), mult. (289->138), div. (0->0), fcn. (220->8), ass. (0->51)
t200 = sin(qJ(1));
t220 = t200 * pkin(1);
t198 = cos(pkin(7));
t219 = t198 * pkin(3);
t195 = pkin(7) + qJ(4);
t191 = sin(t195);
t218 = Icges(5,4) * t191;
t193 = cos(t195);
t217 = Icges(5,4) * t193;
t201 = cos(qJ(1));
t190 = qJD(1) * t201 * pkin(1);
t196 = qJ(1) + pkin(6);
t192 = sin(t196);
t194 = cos(t196);
t215 = qJD(1) * (t194 * pkin(2) + t192 * qJ(3)) + t190;
t214 = qJD(4) * t192;
t213 = qJD(4) * t194;
t212 = -t192 * pkin(2) + t194 * qJ(3) - t220;
t197 = sin(pkin(7));
t211 = rSges(4,1) * t198 - rSges(4,2) * t197;
t210 = rSges(5,1) * t193 - rSges(5,2) * t191;
t209 = Icges(5,1) * t193 - t218;
t208 = -Icges(5,2) * t191 + t217;
t207 = Icges(5,5) * t193 - Icges(5,6) * t191;
t172 = -Icges(5,6) * t194 + t192 * t208;
t174 = -Icges(5,5) * t194 + t192 * t209;
t206 = t172 * t191 - t174 * t193;
t173 = Icges(5,6) * t192 + t194 * t208;
t175 = Icges(5,5) * t192 + t194 * t209;
t205 = -t173 * t191 + t175 * t193;
t182 = Icges(5,2) * t193 + t218;
t183 = Icges(5,1) * t191 + t217;
t204 = -t182 * t191 + t183 * t193;
t202 = qJD(2) ^ 2;
t188 = qJD(3) * t192;
t187 = rSges(2,1) * t201 - rSges(2,2) * t200;
t186 = rSges(2,1) * t200 + rSges(2,2) * t201;
t184 = rSges(5,1) * t191 + rSges(5,2) * t193;
t181 = Icges(5,5) * t191 + Icges(5,6) * t193;
t179 = t190 + qJD(1) * (rSges(3,1) * t194 - rSges(3,2) * t192);
t178 = (-rSges(3,1) * t192 - rSges(3,2) * t194 - t220) * qJD(1);
t177 = t192 * rSges(5,3) + t194 * t210;
t176 = -t194 * rSges(5,3) + t192 * t210;
t171 = Icges(5,3) * t192 + t194 * t207;
t170 = -Icges(5,3) * t194 + t192 * t207;
t169 = qJD(1) * t192 * rSges(4,3) + (qJD(1) * t211 - qJD(3)) * t194 + t215;
t168 = t188 + (t194 * rSges(4,3) - t192 * t211 + t212) * qJD(1);
t167 = qJD(2) + (t176 * t192 + t177 * t194) * qJD(4);
t166 = -t184 * t214 - qJD(3) * t194 + (pkin(5) * t192 + t194 * t219 + t177) * qJD(1) + t215;
t165 = -t184 * t213 + t188 + (pkin(5) * t194 - t192 * t219 - t176 + t212) * qJD(1);
t1 = m(3) * (t178 ^ 2 + t179 ^ 2 + t202) / 0.2e1 + m(4) * (t168 ^ 2 + t169 ^ 2 + t202) / 0.2e1 + m(5) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + ((t192 * t181 + t194 * t204) * qJD(1) + (t192 ^ 2 * t171 + (t206 * t194 + (-t170 + t205) * t192) * t194) * qJD(4)) * t214 / 0.2e1 - ((-t194 * t181 + t192 * t204) * qJD(1) + (t194 ^ 2 * t170 + (t205 * t192 + (-t171 + t206) * t194) * t192) * qJD(4)) * t213 / 0.2e1 + qJD(1) * ((t193 * t182 + t191 * t183) * qJD(1) + ((t193 * t173 + t191 * t175) * t192 - (t172 * t193 + t174 * t191) * t194) * qJD(4)) / 0.2e1 + (m(2) * (t186 ^ 2 + t187 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t198 ^ 2 + (Icges(4,1) * t197 + 0.2e1 * Icges(4,4) * t198) * t197) * qJD(1) ^ 2 / 0.2e1;
T = t1;
