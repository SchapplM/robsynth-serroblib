% Calculate kinetic energy for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:48
% EndTime: 2019-12-31 16:38:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (229->77), mult. (269->127), div. (0->0), fcn. (200->6), ass. (0->45)
t183 = sin(qJ(1));
t200 = t183 * pkin(1);
t182 = sin(qJ(4));
t199 = Icges(5,4) * t182;
t184 = cos(qJ(4));
t198 = Icges(5,4) * t184;
t185 = cos(qJ(1));
t178 = qJD(1) * t185 * pkin(1);
t181 = qJ(1) + pkin(6);
t179 = sin(t181);
t180 = cos(t181);
t197 = qJD(1) * (t180 * pkin(2) + t179 * qJ(3)) + t178;
t196 = qJD(4) * t179;
t195 = -t179 * pkin(2) + t180 * qJ(3) - t200;
t194 = rSges(5,1) * t182 + rSges(5,2) * t184;
t193 = Icges(5,1) * t182 + t198;
t192 = Icges(5,2) * t184 + t199;
t191 = Icges(5,5) * t182 + Icges(5,6) * t184;
t161 = Icges(5,6) * t180 + t192 * t179;
t163 = Icges(5,5) * t180 + t193 * t179;
t190 = -t161 * t184 - t163 * t182;
t162 = Icges(5,6) * t179 - t192 * t180;
t164 = Icges(5,5) * t179 - t193 * t180;
t189 = t162 * t184 + t164 * t182;
t172 = -Icges(5,2) * t182 + t198;
t173 = Icges(5,1) * t184 - t199;
t188 = t172 * t184 + t173 * t182;
t186 = qJD(2) ^ 2;
t177 = qJD(3) * t179;
t176 = t185 * rSges(2,1) - t183 * rSges(2,2);
t175 = t184 * rSges(5,1) - t182 * rSges(5,2);
t174 = t183 * rSges(2,1) + t185 * rSges(2,2);
t171 = Icges(5,5) * t184 - Icges(5,6) * t182;
t168 = t178 + qJD(1) * (t180 * rSges(3,1) - t179 * rSges(3,2));
t167 = (-t179 * rSges(3,1) - t180 * rSges(3,2) - t200) * qJD(1);
t166 = t179 * rSges(5,3) - t194 * t180;
t165 = t180 * rSges(5,3) + t194 * t179;
t160 = Icges(5,3) * t179 - t191 * t180;
t159 = Icges(5,3) * t180 + t191 * t179;
t158 = -qJD(3) * t180 + qJD(1) * (-t180 * rSges(4,2) + t179 * rSges(4,3)) + t197;
t157 = t177 + (t179 * rSges(4,2) + t180 * rSges(4,3) + t195) * qJD(1);
t156 = qJD(2) + (-t165 * t179 + t166 * t180) * qJD(4);
t155 = qJD(1) * t165 + (pkin(5) * qJD(1) - qJD(4) * t175 - qJD(3)) * t180 + t197;
t154 = t175 * t196 + t177 + (-pkin(5) * t179 - t166 + t195) * qJD(1);
t1 = m(3) * (t167 ^ 2 + t168 ^ 2 + t186) / 0.2e1 + m(4) * (t157 ^ 2 + t158 ^ 2 + t186) / 0.2e1 + m(5) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + qJD(4) * t180 * ((t180 * t171 + t188 * t179) * qJD(1) + (t180 ^ 2 * t159 + (t189 * t179 + (t160 - t190) * t180) * t179) * qJD(4)) / 0.2e1 + ((t179 * t171 - t188 * t180) * qJD(1) + (t179 ^ 2 * t160 + (t190 * t180 + (t159 - t189) * t179) * t180) * qJD(4)) * t196 / 0.2e1 + qJD(1) * ((-t182 * t172 + t184 * t173) * qJD(1) + ((-t182 * t161 + t184 * t163) * t180 + (-t182 * t162 + t184 * t164) * t179) * qJD(4)) / 0.2e1 + (m(2) * (t174 ^ 2 + t176 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
