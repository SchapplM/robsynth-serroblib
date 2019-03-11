% Calculate kinetic energy for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:25
% EndTime: 2019-03-08 18:27:25
% DurationCPUTime: 0.05s
% Computational Cost: add. (79->42), mult. (97->64), div. (0->0), fcn. (48->6), ass. (0->25)
t163 = sin(qJ(1));
t170 = t163 * pkin(1);
t165 = cos(qJ(1));
t157 = qJD(1) * t165 * pkin(1);
t161 = qJ(1) + pkin(6);
t158 = sin(t161);
t159 = cos(t161);
t169 = qJD(1) * (t159 * pkin(2) + t158 * qJ(3)) + t157;
t168 = -t158 * pkin(2) + t159 * qJ(3) - t170;
t166 = qJD(2) ^ 2;
t164 = cos(qJ(4));
t162 = sin(qJ(4));
t160 = qJD(1) - qJD(4);
t156 = qJD(3) * t158;
t155 = t165 * rSges(2,1) - t163 * rSges(2,2);
t154 = t163 * rSges(2,1) + t165 * rSges(2,2);
t151 = t158 * t164 - t159 * t162;
t150 = -t158 * t162 - t159 * t164;
t149 = t157 + qJD(1) * (t159 * rSges(3,1) - t158 * rSges(3,2));
t148 = (-t158 * rSges(3,1) - t159 * rSges(3,2) - t170) * qJD(1);
t147 = -qJD(3) * t159 + qJD(1) * (t159 * rSges(4,1) + t158 * rSges(4,3)) + t169;
t146 = t156 + (-t158 * rSges(4,1) + t159 * rSges(4,3) + t168) * qJD(1);
t145 = t160 * (-t150 * rSges(5,1) + t151 * rSges(5,2)) + (qJD(1) * pkin(3) - qJD(3)) * t159 + t169;
t144 = t156 - t160 * (t151 * rSges(5,1) + t150 * rSges(5,2)) + (-t158 * pkin(3) + t168) * qJD(1);
t1 = m(3) * (t148 ^ 2 + t149 ^ 2 + t166) / 0.2e1 + m(4) * (t146 ^ 2 + t147 ^ 2 + t166) / 0.2e1 + m(5) * (t144 ^ 2 + t145 ^ 2 + t166) / 0.2e1 + t160 ^ 2 * Icges(5,3) / 0.2e1 + (m(2) * (t154 ^ 2 + t155 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
