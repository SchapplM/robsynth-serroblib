% Calculate kinetic energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:47
% EndTime: 2019-12-31 16:20:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (320->75), mult. (267->130), div. (0->0), fcn. (210->6), ass. (0->46)
t196 = cos(pkin(7));
t213 = t196 * pkin(3);
t193 = pkin(7) + qJ(4);
t189 = sin(t193);
t212 = Icges(5,4) * t189;
t191 = cos(t193);
t211 = Icges(5,4) * t191;
t194 = pkin(6) + qJ(2);
t190 = sin(t194);
t209 = qJD(4) * t190;
t192 = cos(t194);
t208 = qJD(4) * t192;
t195 = sin(pkin(7));
t207 = rSges(4,1) * t196 - rSges(4,2) * t195;
t206 = rSges(5,1) * t191 - rSges(5,2) * t189;
t205 = Icges(5,1) * t191 - t212;
t204 = -Icges(5,2) * t189 + t211;
t203 = Icges(5,5) * t191 - Icges(5,6) * t189;
t173 = -Icges(5,6) * t192 + t190 * t204;
t175 = -Icges(5,5) * t192 + t190 * t205;
t202 = t173 * t189 - t175 * t191;
t174 = Icges(5,6) * t190 + t192 * t204;
t176 = Icges(5,5) * t190 + t192 * t205;
t201 = -t174 * t189 + t176 * t191;
t181 = Icges(5,2) * t191 + t212;
t182 = Icges(5,1) * t189 + t211;
t200 = -t181 * t189 + t182 * t191;
t199 = qJD(1) ^ 2;
t198 = qJD(2) ^ 2;
t187 = qJD(3) * t190;
t186 = t192 * rSges(3,1) - t190 * rSges(3,2);
t185 = t190 * rSges(3,1) + t192 * rSges(3,2);
t184 = t190 * pkin(2) - t192 * qJ(3);
t183 = t189 * rSges(5,1) + t191 * rSges(5,2);
t180 = Icges(5,5) * t189 + Icges(5,6) * t191;
t179 = qJD(2) * (t192 * pkin(2) + t190 * qJ(3));
t178 = t190 * rSges(5,3) + t192 * t206;
t177 = -t192 * rSges(5,3) + t190 * t206;
t172 = Icges(5,3) * t190 + t192 * t203;
t171 = -Icges(5,3) * t192 + t190 * t203;
t170 = qJD(2) * t190 * rSges(4,3) + t179 + (qJD(2) * t207 - qJD(3)) * t192;
t169 = t187 + (t192 * rSges(4,3) - t190 * t207 - t184) * qJD(2);
t168 = qJD(1) + (t177 * t190 + t178 * t192) * qJD(4);
t167 = -t183 * t209 - qJD(3) * t192 + t179 + (pkin(5) * t190 + t192 * t213 + t178) * qJD(2);
t166 = -t183 * t208 + t187 + (pkin(5) * t192 - t190 * t213 - t177 - t184) * qJD(2);
t1 = m(2) * t199 / 0.2e1 + m(3) * (t199 + (t185 ^ 2 + t186 ^ 2) * t198) / 0.2e1 + m(4) * (t169 ^ 2 + t170 ^ 2 + t199) / 0.2e1 + m(5) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + ((t190 * t180 + t192 * t200) * qJD(2) + (t190 ^ 2 * t172 + (t202 * t192 + (-t171 + t201) * t190) * t192) * qJD(4)) * t209 / 0.2e1 - ((-t192 * t180 + t190 * t200) * qJD(2) + (t192 ^ 2 * t171 + (t201 * t190 + (-t172 + t202) * t192) * t190) * qJD(4)) * t208 / 0.2e1 + qJD(2) * ((t191 * t181 + t189 * t182) * qJD(2) + ((t191 * t174 + t189 * t176) * t190 - (t191 * t173 + t189 * t175) * t192) * qJD(4)) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t196 ^ 2 + (Icges(4,1) * t195 + 0.2e1 * Icges(4,4) * t196) * t195) * t198 / 0.2e1;
T = t1;
