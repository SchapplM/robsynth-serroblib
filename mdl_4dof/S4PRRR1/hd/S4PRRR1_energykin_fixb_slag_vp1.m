% Calculate kinetic energy for
% S4PRRR1
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:04
% EndTime: 2019-03-08 18:25:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (67->29), mult. (59->47), div. (0->0), fcn. (18->6), ass. (0->24)
t154 = qJD(2) + qJD(3);
t159 = pkin(3) * t154;
t158 = pkin(2) * qJD(2);
t153 = pkin(7) + qJ(2);
t149 = sin(t153);
t157 = t149 * t158;
t152 = qJ(3) + t153;
t156 = qJD(1) ^ 2;
t155 = qJD(2) ^ 2;
t151 = qJD(4) + t154;
t150 = cos(t153);
t148 = qJ(4) + t152;
t147 = cos(t152);
t146 = sin(t152);
t145 = cos(t148);
t144 = sin(t148);
t143 = t150 * t158;
t142 = rSges(3,1) * t150 - rSges(3,2) * t149;
t141 = rSges(3,1) * t149 + rSges(3,2) * t150;
t140 = t143 + t154 * (rSges(4,1) * t147 - rSges(4,2) * t146);
t139 = -t157 - t154 * (rSges(4,1) * t146 + rSges(4,2) * t147);
t138 = t143 + t147 * t159 + t151 * (rSges(5,1) * t145 - rSges(5,2) * t144);
t137 = -t157 - t146 * t159 - t151 * (rSges(5,1) * t144 + rSges(5,2) * t145);
t1 = m(2) * t156 / 0.2e1 + m(3) * (t156 + (t141 ^ 2 + t142 ^ 2) * t155) / 0.2e1 + t155 * Icges(3,3) / 0.2e1 + m(4) * (t139 ^ 2 + t140 ^ 2 + t156) / 0.2e1 + t154 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t156) / 0.2e1 + t151 ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
