% Calculate kinetic energy for
% S4RRPR1
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
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:54
% EndTime: 2019-03-08 18:34:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->37), mult. (80->59), div. (0->0), fcn. (28->8), ass. (0->28)
t159 = qJ(1) + qJ(2);
t155 = sin(t159);
t167 = pkin(2) * t155;
t166 = pkin(1) * qJD(1);
t161 = cos(qJ(1));
t152 = t161 * t166;
t156 = cos(t159);
t158 = qJD(1) + qJD(2);
t165 = t158 * pkin(2) * t156 + t152;
t160 = sin(qJ(1));
t164 = t160 * t166;
t154 = pkin(7) + t159;
t162 = qJD(3) ^ 2;
t153 = qJD(4) + t158;
t151 = qJ(4) + t154;
t150 = cos(t154);
t149 = sin(t154);
t148 = cos(t151);
t147 = sin(t151);
t146 = rSges(2,1) * t161 - rSges(2,2) * t160;
t145 = rSges(2,1) * t160 + rSges(2,2) * t161;
t143 = t152 + t158 * (rSges(3,1) * t156 - rSges(3,2) * t155);
t142 = -t164 - t158 * (rSges(3,1) * t155 + rSges(3,2) * t156);
t141 = t158 * (rSges(4,1) * t150 - rSges(4,2) * t149) + t165;
t140 = -t164 + (-rSges(4,1) * t149 - rSges(4,2) * t150 - t167) * t158;
t139 = t158 * pkin(3) * t150 + t153 * (rSges(5,1) * t148 - rSges(5,2) * t147) + t165;
t138 = -t164 - t153 * (rSges(5,1) * t147 + rSges(5,2) * t148) + (-pkin(3) * t149 - t167) * t158;
t1 = m(3) * (t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t162) / 0.2e1 + m(5) * (t138 ^ 2 + t139 ^ 2 + t162) / 0.2e1 + t153 ^ 2 * Icges(5,3) / 0.2e1 + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t158 ^ 2 + (m(2) * (t145 ^ 2 + t146 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
