% Calculate kinetic energy for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:27
% EndTime: 2019-12-31 17:03:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (244->75), mult. (267->130), div. (0->0), fcn. (200->6), ass. (0->46)
t202 = pkin(1) * qJD(1);
t184 = sin(qJ(4));
t201 = Icges(5,4) * t184;
t186 = cos(qJ(4));
t200 = Icges(5,4) * t186;
t187 = cos(qJ(1));
t178 = t187 * t202;
t183 = qJ(1) + qJ(2);
t179 = sin(t183);
t180 = cos(t183);
t182 = qJD(1) + qJD(2);
t199 = t182 * (t180 * pkin(2) + t179 * qJ(3)) + t178;
t198 = qJD(4) * t179;
t185 = sin(qJ(1));
t197 = t185 * t202;
t196 = qJD(3) * t179 - t197;
t195 = rSges(5,1) * t184 + rSges(5,2) * t186;
t194 = Icges(5,1) * t184 + t200;
t193 = Icges(5,2) * t186 + t201;
t192 = Icges(5,5) * t184 + Icges(5,6) * t186;
t161 = Icges(5,6) * t180 + t193 * t179;
t163 = Icges(5,5) * t180 + t194 * t179;
t191 = -t161 * t186 - t163 * t184;
t162 = Icges(5,6) * t179 - t193 * t180;
t164 = Icges(5,5) * t179 - t194 * t180;
t190 = t162 * t186 + t164 * t184;
t172 = -Icges(5,2) * t184 + t200;
t173 = Icges(5,1) * t186 - t201;
t189 = t172 * t186 + t173 * t184;
t176 = t187 * rSges(2,1) - t185 * rSges(2,2);
t175 = t186 * rSges(5,1) - t184 * rSges(5,2);
t174 = t185 * rSges(2,1) + t187 * rSges(2,2);
t171 = Icges(5,5) * t186 - Icges(5,6) * t184;
t170 = t179 * pkin(2) - t180 * qJ(3);
t168 = t178 + t182 * (t180 * rSges(3,1) - t179 * rSges(3,2));
t167 = -t197 - t182 * (t179 * rSges(3,1) + t180 * rSges(3,2));
t166 = t179 * rSges(5,3) - t195 * t180;
t165 = t180 * rSges(5,3) + t195 * t179;
t160 = Icges(5,3) * t179 - t192 * t180;
t159 = Icges(5,3) * t180 + t192 * t179;
t158 = -qJD(3) * t180 + t182 * (-t180 * rSges(4,2) + t179 * rSges(4,3)) + t199;
t157 = (t179 * rSges(4,2) + t180 * rSges(4,3) - t170) * t182 + t196;
t156 = (-t165 * t179 + t166 * t180) * qJD(4);
t155 = t182 * t165 + (pkin(6) * t182 - qJD(4) * t175 - qJD(3)) * t180 + t199;
t154 = t175 * t198 + (-pkin(6) * t179 - t166 - t170) * t182 + t196;
t1 = m(3) * (t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(4) * (t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(5) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + qJD(4) * t180 * ((t180 * t171 + t189 * t179) * t182 + (t180 ^ 2 * t159 + (t190 * t179 + (t160 - t191) * t180) * t179) * qJD(4)) / 0.2e1 + ((t179 * t171 - t189 * t180) * t182 + (t179 ^ 2 * t160 + (t191 * t180 + (t159 - t190) * t179) * t180) * qJD(4)) * t198 / 0.2e1 + t182 * ((-t184 * t172 + t186 * t173) * t182 + ((-t184 * t161 + t186 * t163) * t180 + (-t184 * t162 + t186 * t164) * t179) * qJD(4)) / 0.2e1 + (Icges(3,3) + Icges(4,1)) * t182 ^ 2 / 0.2e1 + (m(2) * (t174 ^ 2 + t176 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
