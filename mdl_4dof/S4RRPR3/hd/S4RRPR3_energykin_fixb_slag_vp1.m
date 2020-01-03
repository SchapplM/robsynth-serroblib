% Calculate kinetic energy for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:27
% EndTime: 2019-12-31 17:01:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (332->75), mult. (262->128), div. (0->0), fcn. (194->8), ass. (0->49)
t181 = qJ(1) + qJ(2);
t177 = sin(t181);
t201 = pkin(2) * t177;
t200 = pkin(1) * qJD(1);
t182 = sin(qJ(4));
t199 = Icges(5,4) * t182;
t184 = cos(qJ(4));
t198 = Icges(5,4) * t184;
t185 = cos(qJ(1));
t175 = t185 * t200;
t178 = cos(t181);
t180 = qJD(1) + qJD(2);
t197 = t180 * pkin(2) * t178 + t175;
t176 = pkin(7) + t181;
t173 = sin(t176);
t196 = qJD(4) * t173;
t174 = cos(t176);
t195 = qJD(4) * t174;
t183 = sin(qJ(1));
t194 = t183 * t200;
t193 = rSges(5,1) * t184 - rSges(5,2) * t182;
t192 = Icges(5,1) * t184 - t199;
t191 = -Icges(5,2) * t182 + t198;
t190 = Icges(5,5) * t184 - Icges(5,6) * t182;
t158 = -Icges(5,6) * t174 + t191 * t173;
t160 = -Icges(5,5) * t174 + t192 * t173;
t189 = t158 * t182 - t160 * t184;
t159 = Icges(5,6) * t173 + t191 * t174;
t161 = Icges(5,5) * t173 + t192 * t174;
t188 = -t159 * t182 + t161 * t184;
t167 = Icges(5,2) * t184 + t199;
t168 = Icges(5,1) * t182 + t198;
t187 = -t167 * t182 + t168 * t184;
t172 = rSges(2,1) * t185 - rSges(2,2) * t183;
t171 = rSges(2,1) * t183 + rSges(2,2) * t185;
t170 = rSges(5,1) * t182 + rSges(5,2) * t184;
t166 = Icges(5,5) * t182 + Icges(5,6) * t184;
t165 = t175 + t180 * (rSges(3,1) * t178 - rSges(3,2) * t177);
t164 = -t194 - t180 * (rSges(3,1) * t177 + rSges(3,2) * t178);
t163 = t173 * rSges(5,3) + t193 * t174;
t162 = -t174 * rSges(5,3) + t193 * t173;
t157 = Icges(5,3) * t173 + t190 * t174;
t156 = -Icges(5,3) * t174 + t190 * t173;
t155 = t180 * (rSges(4,1) * t174 - rSges(4,2) * t173) + t197;
t154 = -t194 + (-rSges(4,1) * t173 - rSges(4,2) * t174 - t201) * t180;
t153 = qJD(3) + (t162 * t173 + t163 * t174) * qJD(4);
t152 = -t170 * t196 + (pkin(3) * t174 + pkin(6) * t173 + t163) * t180 + t197;
t151 = -t194 - t170 * t195 + (-pkin(3) * t173 + pkin(6) * t174 - t162 - t201) * t180;
t1 = m(3) * (t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(5) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + ((t173 * t166 + t187 * t174) * t180 + (t173 ^ 2 * t157 + (t189 * t174 + (-t156 + t188) * t173) * t174) * qJD(4)) * t196 / 0.2e1 - ((-t174 * t166 + t187 * t173) * t180 + (t174 ^ 2 * t156 + (t188 * t173 + (-t157 + t189) * t174) * t173) * qJD(4)) * t195 / 0.2e1 + t180 * ((t184 * t167 + t182 * t168) * t180 + ((t159 * t184 + t161 * t182) * t173 - (t158 * t184 + t160 * t182) * t174) * qJD(4)) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t180 ^ 2 / 0.2e1 + (m(2) * (t171 ^ 2 + t172 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
