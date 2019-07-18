% Calculate kinetic energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:29
% EndTime: 2019-07-18 18:16:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (91->40), mult. (94->66), div. (0->0), fcn. (48->6), ass. (0->26)
t169 = pkin(1) * qJD(1);
t164 = cos(qJ(1));
t154 = t164 * t169;
t160 = qJ(1) + qJ(2);
t156 = sin(t160);
t157 = cos(t160);
t159 = qJD(1) + qJD(2);
t168 = t159 * (t157 * pkin(2) + t156 * qJ(3)) + t154;
t162 = sin(qJ(1));
t167 = t162 * t169;
t166 = qJD(3) * t156 - t167;
t163 = cos(qJ(4));
t161 = sin(qJ(4));
t155 = -qJD(4) + t159;
t152 = t164 * rSges(2,1) - t162 * rSges(2,2);
t151 = t162 * rSges(2,1) + t164 * rSges(2,2);
t150 = t156 * pkin(2) - t157 * qJ(3);
t149 = t156 * t163 - t157 * t161;
t148 = -t156 * t161 - t157 * t163;
t146 = t154 + t159 * (t157 * rSges(3,1) - t156 * rSges(3,2));
t145 = -t167 - t159 * (t156 * rSges(3,1) + t157 * rSges(3,2));
t144 = -qJD(3) * t157 + t159 * (t157 * rSges(4,1) + t156 * rSges(4,3)) + t168;
t143 = (-t156 * rSges(4,1) + t157 * rSges(4,3) - t150) * t159 + t166;
t142 = t155 * (-t148 * rSges(5,1) + t149 * rSges(5,2)) + (t159 * pkin(3) - qJD(3)) * t157 + t168;
t141 = -t155 * (t149 * rSges(5,1) + t148 * rSges(5,2)) + (-t156 * pkin(3) - t150) * t159 + t166;
t1 = m(3) * (t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(4) * (t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2) / 0.2e1 + t155 ^ 2 * Icges(5,3) / 0.2e1 + (Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t159 ^ 2 + (m(2) * (t151 ^ 2 + t152 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
