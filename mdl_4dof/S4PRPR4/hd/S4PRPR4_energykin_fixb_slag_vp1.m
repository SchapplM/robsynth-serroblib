% Calculate kinetic energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:55
% EndTime: 2019-12-31 16:21:55
% DurationCPUTime: 0.22s
% Computational Cost: add. (219->69), mult. (247->119), div. (0->0), fcn. (190->4), ass. (0->40)
t182 = sin(qJ(4));
t195 = Icges(5,4) * t182;
t183 = cos(qJ(4));
t194 = Icges(5,4) * t183;
t181 = pkin(6) + qJ(2);
t179 = sin(t181);
t193 = qJD(4) * t179;
t192 = rSges(5,1) * t182 + rSges(5,2) * t183;
t191 = Icges(5,1) * t182 + t194;
t190 = Icges(5,2) * t183 + t195;
t189 = Icges(5,5) * t182 + Icges(5,6) * t183;
t180 = cos(t181);
t164 = Icges(5,6) * t180 + t179 * t190;
t166 = Icges(5,5) * t180 + t179 * t191;
t188 = -t164 * t183 - t166 * t182;
t165 = Icges(5,6) * t179 - t180 * t190;
t167 = Icges(5,5) * t179 - t180 * t191;
t187 = t165 * t183 + t167 * t182;
t175 = -Icges(5,2) * t182 + t194;
t176 = Icges(5,1) * t183 - t195;
t186 = t175 * t183 + t176 * t182;
t185 = qJD(1) ^ 2;
t184 = qJD(2) ^ 2;
t178 = qJD(3) * t179;
t177 = t183 * rSges(5,1) - t182 * rSges(5,2);
t174 = Icges(5,5) * t183 - Icges(5,6) * t182;
t173 = t180 * rSges(3,1) - t179 * rSges(3,2);
t172 = t179 * rSges(3,1) + t180 * rSges(3,2);
t171 = t179 * pkin(2) - t180 * qJ(3);
t170 = qJD(2) * (t180 * pkin(2) + t179 * qJ(3));
t169 = t179 * rSges(5,3) - t180 * t192;
t168 = t180 * rSges(5,3) + t179 * t192;
t163 = Icges(5,3) * t179 - t180 * t189;
t162 = Icges(5,3) * t180 + t179 * t189;
t161 = t170 - qJD(3) * t180 + qJD(2) * (-t180 * rSges(4,2) + t179 * rSges(4,3));
t160 = t178 + (t179 * rSges(4,2) + t180 * rSges(4,3) - t171) * qJD(2);
t159 = qJD(1) + (-t168 * t179 + t169 * t180) * qJD(4);
t158 = qJD(2) * t168 + t170 + (pkin(5) * qJD(2) - qJD(4) * t177 - qJD(3)) * t180;
t157 = t177 * t193 + t178 + (-pkin(5) * t179 - t169 - t171) * qJD(2);
t1 = m(2) * t185 / 0.2e1 + m(3) * (t185 + (t172 ^ 2 + t173 ^ 2) * t184) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t185) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + qJD(4) * t180 * ((t180 * t174 + t179 * t186) * qJD(2) + (t180 ^ 2 * t162 + (t187 * t179 + (t163 - t188) * t180) * t179) * qJD(4)) / 0.2e1 + ((t179 * t174 - t180 * t186) * qJD(2) + (t179 ^ 2 * t163 + (t188 * t180 + (t162 - t187) * t179) * t180) * qJD(4)) * t193 / 0.2e1 + qJD(2) * ((-t182 * t175 + t183 * t176) * qJD(2) + ((-t182 * t164 + t183 * t166) * t180 + (-t182 * t165 + t183 * t167) * t179) * qJD(4)) / 0.2e1 + (Icges(3,3) + Icges(4,1)) * t184 / 0.2e1;
T = t1;
