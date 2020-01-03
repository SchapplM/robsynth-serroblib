% Calculate kinetic energy for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:33
% EndTime: 2019-12-31 16:31:34
% DurationCPUTime: 0.17s
% Computational Cost: add. (316->67), mult. (241->119), div. (0->0), fcn. (184->6), ass. (0->45)
t195 = pkin(2) * qJD(2);
t179 = sin(qJ(4));
t194 = Icges(5,4) * t179;
t180 = cos(qJ(4));
t193 = Icges(5,4) * t180;
t177 = pkin(7) + qJ(2);
t176 = qJ(3) + t177;
t172 = sin(t176);
t192 = qJD(4) * t172;
t173 = cos(t176);
t191 = qJD(4) * t173;
t174 = sin(t177);
t190 = t174 * t195;
t189 = rSges(5,1) * t180 - rSges(5,2) * t179;
t188 = Icges(5,1) * t180 - t194;
t187 = -Icges(5,2) * t179 + t193;
t186 = Icges(5,5) * t180 - Icges(5,6) * t179;
t157 = -Icges(5,6) * t173 + t187 * t172;
t159 = -Icges(5,5) * t173 + t188 * t172;
t185 = t157 * t179 - t159 * t180;
t158 = Icges(5,6) * t172 + t187 * t173;
t160 = Icges(5,5) * t172 + t188 * t173;
t184 = -t158 * t179 + t160 * t180;
t168 = Icges(5,2) * t180 + t194;
t169 = Icges(5,1) * t179 + t193;
t183 = -t168 * t179 + t169 * t180;
t182 = qJD(1) ^ 2;
t181 = qJD(2) ^ 2;
t178 = qJD(2) + qJD(3);
t175 = cos(t177);
t171 = t175 * t195;
t170 = t179 * rSges(5,1) + t180 * rSges(5,2);
t167 = Icges(5,5) * t179 + Icges(5,6) * t180;
t166 = t175 * rSges(3,1) - t174 * rSges(3,2);
t165 = t174 * rSges(3,1) + t175 * rSges(3,2);
t164 = t171 + t178 * (t173 * rSges(4,1) - t172 * rSges(4,2));
t163 = -t190 - t178 * (t172 * rSges(4,1) + t173 * rSges(4,2));
t162 = t172 * rSges(5,3) + t189 * t173;
t161 = -t173 * rSges(5,3) + t189 * t172;
t156 = Icges(5,3) * t172 + t186 * t173;
t155 = -Icges(5,3) * t173 + t186 * t172;
t154 = qJD(1) + (t161 * t172 + t162 * t173) * qJD(4);
t153 = -t170 * t192 + t171 + (t173 * pkin(3) + t172 * pkin(6) + t162) * t178;
t152 = -t190 - t170 * t191 + (-t172 * pkin(3) + t173 * pkin(6) - t161) * t178;
t1 = m(2) * t182 / 0.2e1 + m(3) * (t182 + (t165 ^ 2 + t166 ^ 2) * t181) / 0.2e1 + t181 * Icges(3,3) / 0.2e1 + m(4) * (t163 ^ 2 + t164 ^ 2 + t182) / 0.2e1 + t178 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + ((t172 * t167 + t183 * t173) * t178 + (t172 ^ 2 * t156 + (t185 * t173 + (-t155 + t184) * t172) * t173) * qJD(4)) * t192 / 0.2e1 - ((-t173 * t167 + t183 * t172) * t178 + (t173 ^ 2 * t155 + (t184 * t172 + (-t156 + t185) * t173) * t172) * qJD(4)) * t191 / 0.2e1 + t178 * ((t180 * t168 + t179 * t169) * t178 + ((t180 * t158 + t179 * t160) * t172 - (t180 * t157 + t179 * t159) * t173) * qJD(4)) / 0.2e1;
T = t1;
