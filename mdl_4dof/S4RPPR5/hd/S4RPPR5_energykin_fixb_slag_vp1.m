% Calculate kinetic energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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

function T = S4RPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:42
% EndTime: 2019-12-31 16:39:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (220->61), mult. (470->107), div. (0->0), fcn. (500->6), ass. (0->35)
t209 = sin(pkin(6));
t210 = cos(pkin(6));
t211 = sin(qJ(1));
t212 = cos(qJ(1));
t178 = -t209 * t211 - t210 * t212;
t219 = t178 ^ 2;
t179 = t209 * t212 - t210 * t211;
t218 = t179 ^ 2;
t193 = sin(qJ(4));
t194 = cos(qJ(4));
t217 = -Icges(5,5) * t193 - Icges(5,6) * t194;
t215 = t179 * t178;
t214 = t217 * qJD(1);
t206 = qJD(4) * (-t193 * rSges(5,1) - t194 * rSges(5,2));
t185 = pkin(1) * t211 - qJ(2) * t212;
t205 = -pkin(2) * t211 - t185;
t204 = -qJD(2) * t212 + qJD(1) * (pkin(1) * t212 + qJ(2) * t211);
t203 = -rSges(5,1) * t194 + rSges(5,2) * t193;
t200 = -Icges(5,5) * t194 + Icges(5,6) * t193;
t196 = qJD(1) * t212 * pkin(2) + t204;
t192 = qJD(2) * t211;
t187 = rSges(2,1) * t212 - rSges(2,2) * t211;
t186 = rSges(2,1) * t211 + rSges(2,2) * t212;
t177 = qJD(1) * (rSges(3,1) * t212 + rSges(3,3) * t211) + t204;
t176 = t192 + (-rSges(3,1) * t211 + rSges(3,3) * t212 - t185) * qJD(1);
t175 = t179 * rSges(5,3) + t178 * t203;
t174 = -t178 * rSges(5,3) + t179 * t203;
t169 = Icges(5,3) * t179 + t178 * t200;
t168 = -Icges(5,3) * t178 + t179 * t200;
t167 = qJD(1) * (-t178 * rSges(4,1) - t179 * rSges(4,2)) + t196;
t166 = t192 + (t179 * rSges(4,1) - t178 * rSges(4,2) + t205) * qJD(1);
t165 = -t179 * t206 + (-t178 * pkin(3) + t179 * pkin(5) + t175) * qJD(1) + t196;
t164 = -t178 * t206 + t192 + (t179 * pkin(3) + t178 * pkin(5) - t174 + t205) * qJD(1);
t163 = -qJD(3) + (t174 * t179 + t175 * t178) * qJD(4);
t1 = m(3) * (t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(5) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + qJD(4) * t179 * (t179 * t214 + (-t168 * t215 + t218 * t169) * qJD(4)) / 0.2e1 - qJD(4) * t178 * (-t178 * t214 + (t219 * t168 - t169 * t215) * qJD(4)) / 0.2e1 + qJD(1) * ((t194 ^ 2 * Icges(5,2) + (Icges(5,1) * t193 + 0.2e1 * Icges(5,4) * t194) * t193) * qJD(1) + (t218 + t219) * t217 * qJD(4)) / 0.2e1 + (m(2) * (t186 ^ 2 + t187 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
