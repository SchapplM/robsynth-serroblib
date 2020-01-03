% Calculate kinetic energy for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:07
% DurationCPUTime: 0.24s
% Computational Cost: add. (326->75), mult. (263->128), div. (0->0), fcn. (194->8), ass. (0->49)
t183 = sin(qJ(1));
t201 = pkin(1) * t183;
t182 = sin(qJ(4));
t200 = Icges(5,4) * t182;
t184 = cos(qJ(4));
t199 = Icges(5,4) * t184;
t185 = cos(qJ(1));
t176 = qJD(1) * t185 * pkin(1);
t181 = qJ(1) + pkin(7);
t178 = cos(t181);
t198 = qJD(1) * pkin(2) * t178 + t176;
t179 = qJ(3) + t181;
t174 = sin(t179);
t197 = qJD(4) * t174;
t175 = cos(t179);
t196 = qJD(4) * t175;
t195 = rSges(5,1) * t184 - rSges(5,2) * t182;
t194 = Icges(5,1) * t184 - t200;
t193 = -Icges(5,2) * t182 + t199;
t192 = Icges(5,5) * t184 - Icges(5,6) * t182;
t159 = -Icges(5,6) * t175 + t193 * t174;
t161 = -Icges(5,5) * t175 + t194 * t174;
t191 = t159 * t182 - t161 * t184;
t160 = Icges(5,6) * t174 + t193 * t175;
t162 = Icges(5,5) * t174 + t194 * t175;
t190 = -t160 * t182 + t162 * t184;
t168 = Icges(5,2) * t184 + t200;
t169 = Icges(5,1) * t182 + t199;
t189 = -t168 * t182 + t169 * t184;
t177 = sin(t181);
t188 = (-pkin(2) * t177 - t201) * qJD(1);
t186 = qJD(2) ^ 2;
t180 = qJD(1) + qJD(3);
t172 = rSges(2,1) * t185 - rSges(2,2) * t183;
t171 = rSges(2,1) * t183 + rSges(2,2) * t185;
t170 = rSges(5,1) * t182 + rSges(5,2) * t184;
t167 = Icges(5,5) * t182 + Icges(5,6) * t184;
t166 = t176 + qJD(1) * (rSges(3,1) * t178 - rSges(3,2) * t177);
t165 = (-rSges(3,1) * t177 - rSges(3,2) * t178 - t201) * qJD(1);
t164 = t174 * rSges(5,3) + t195 * t175;
t163 = -t175 * rSges(5,3) + t195 * t174;
t158 = Icges(5,3) * t174 + t192 * t175;
t157 = -Icges(5,3) * t175 + t192 * t174;
t156 = t180 * (rSges(4,1) * t175 - rSges(4,2) * t174) + t198;
t155 = -t180 * (rSges(4,1) * t174 + rSges(4,2) * t175) + t188;
t154 = qJD(2) + (t163 * t174 + t164 * t175) * qJD(4);
t153 = -t170 * t197 + (pkin(3) * t175 + pkin(6) * t174 + t164) * t180 + t198;
t152 = -t170 * t196 + (-pkin(3) * t174 + pkin(6) * t175 - t163) * t180 + t188;
t1 = m(3) * (t165 ^ 2 + t166 ^ 2 + t186) / 0.2e1 + m(4) * (t155 ^ 2 + t156 ^ 2 + t186) / 0.2e1 + t180 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + ((t174 * t167 + t189 * t175) * t180 + (t174 ^ 2 * t158 + (t191 * t175 + (-t157 + t190) * t174) * t175) * qJD(4)) * t197 / 0.2e1 - ((-t175 * t167 + t189 * t174) * t180 + (t175 ^ 2 * t157 + (t190 * t174 + (-t158 + t191) * t175) * t174) * qJD(4)) * t196 / 0.2e1 + t180 * ((t184 * t168 + t182 * t169) * t180 + ((t160 * t184 + t162 * t182) * t174 - (t159 * t184 + t161 * t182) * t175) * qJD(4)) / 0.2e1 + (m(2) * (t171 ^ 2 + t172 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
