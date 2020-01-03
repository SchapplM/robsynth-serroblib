% Calculate kinetic energy for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:07
% EndTime: 2019-12-31 17:22:07
% DurationCPUTime: 0.23s
% Computational Cost: add. (341->73), mult. (261->130), div. (0->0), fcn. (194->8), ass. (0->51)
t179 = qJD(1) + qJD(2);
t201 = pkin(2) * t179;
t200 = pkin(1) * qJD(1);
t181 = sin(qJ(4));
t199 = Icges(5,4) * t181;
t183 = cos(qJ(4));
t198 = Icges(5,4) * t183;
t180 = qJ(1) + qJ(2);
t184 = cos(qJ(1));
t174 = t184 * t200;
t177 = cos(t180);
t197 = t177 * t201 + t174;
t178 = qJ(3) + t180;
t172 = sin(t178);
t196 = qJD(4) * t172;
t173 = cos(t178);
t195 = qJD(4) * t173;
t182 = sin(qJ(1));
t194 = t182 * t200;
t193 = rSges(5,1) * t183 - rSges(5,2) * t181;
t192 = Icges(5,1) * t183 - t199;
t191 = -Icges(5,2) * t181 + t198;
t190 = Icges(5,5) * t183 - Icges(5,6) * t181;
t157 = -Icges(5,6) * t173 + t191 * t172;
t159 = -Icges(5,5) * t173 + t192 * t172;
t189 = t157 * t181 - t159 * t183;
t158 = Icges(5,6) * t172 + t191 * t173;
t160 = Icges(5,5) * t172 + t192 * t173;
t188 = -t158 * t181 + t160 * t183;
t166 = Icges(5,2) * t183 + t199;
t167 = Icges(5,1) * t181 + t198;
t187 = -t166 * t181 + t167 * t183;
t176 = sin(t180);
t186 = -t176 * t201 - t194;
t175 = qJD(3) + t179;
t171 = rSges(2,1) * t184 - rSges(2,2) * t182;
t170 = rSges(2,1) * t182 + rSges(2,2) * t184;
t169 = rSges(5,1) * t181 + rSges(5,2) * t183;
t165 = Icges(5,5) * t181 + Icges(5,6) * t183;
t164 = t174 + t179 * (rSges(3,1) * t177 - rSges(3,2) * t176);
t163 = -t194 - t179 * (rSges(3,1) * t176 + rSges(3,2) * t177);
t162 = t172 * rSges(5,3) + t193 * t173;
t161 = -t173 * rSges(5,3) + t193 * t172;
t156 = Icges(5,3) * t172 + t190 * t173;
t155 = -Icges(5,3) * t173 + t190 * t172;
t154 = t175 * (rSges(4,1) * t173 - rSges(4,2) * t172) + t197;
t153 = -t175 * (rSges(4,1) * t172 + rSges(4,2) * t173) + t186;
t152 = (t161 * t172 + t162 * t173) * qJD(4);
t151 = -t169 * t196 + (pkin(3) * t173 + pkin(7) * t172 + t162) * t175 + t197;
t150 = -t169 * t195 + (-pkin(3) * t172 + pkin(7) * t173 - t161) * t175 + t186;
t1 = m(3) * (t163 ^ 2 + t164 ^ 2) / 0.2e1 + t179 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t153 ^ 2 + t154 ^ 2) / 0.2e1 + t175 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + ((t172 * t165 + t187 * t173) * t175 + (t172 ^ 2 * t156 + (t189 * t173 + (-t155 + t188) * t172) * t173) * qJD(4)) * t196 / 0.2e1 - ((-t173 * t165 + t187 * t172) * t175 + (t173 ^ 2 * t155 + (t188 * t172 + (-t156 + t189) * t173) * t172) * qJD(4)) * t195 / 0.2e1 + t175 * ((t183 * t166 + t181 * t167) * t175 + ((t158 * t183 + t160 * t181) * t172 - (t157 * t183 + t159 * t181) * t173) * qJD(4)) / 0.2e1 + (m(2) * (t170 ^ 2 + t171 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
