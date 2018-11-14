% Calculate kinetic energy for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:28
% EndTime: 2018-11-14 13:50:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->37), mult. (81->58), div. (0->0), fcn. (28->8), ass. (0->28)
t158 = sin(qJ(1));
t165 = pkin(1) * t158;
t156 = qJD(1) + qJD(3);
t164 = pkin(3) * t156;
t159 = cos(qJ(1));
t151 = qJD(1) * t159 * pkin(1);
t157 = qJ(1) + pkin(7);
t153 = cos(t157);
t163 = qJD(1) * pkin(2) * t153 + t151;
t155 = qJ(3) + t157;
t152 = sin(t157);
t162 = (-pkin(2) * t152 - t165) * qJD(1);
t160 = qJD(2) ^ 2;
t154 = qJD(4) + t156;
t150 = qJ(4) + t155;
t149 = cos(t155);
t148 = sin(t155);
t147 = cos(t150);
t146 = sin(t150);
t144 = rSges(2,1) * t159 - rSges(2,2) * t158;
t143 = rSges(2,1) * t158 + rSges(2,2) * t159;
t142 = t151 + qJD(1) * (rSges(3,1) * t153 - rSges(3,2) * t152);
t141 = (-rSges(3,1) * t152 - rSges(3,2) * t153 - t165) * qJD(1);
t140 = t156 * (rSges(4,1) * t149 - rSges(4,2) * t148) + t163;
t139 = -t156 * (rSges(4,1) * t148 + rSges(4,2) * t149) + t162;
t138 = t149 * t164 + t154 * (rSges(5,1) * t147 - rSges(5,2) * t146) + t163;
t137 = -t148 * t164 - t154 * (rSges(5,1) * t146 + rSges(5,2) * t147) + t162;
t1 = m(3) * (t141 ^ 2 + t142 ^ 2 + t160) / 0.2e1 + m(4) * (t139 ^ 2 + t140 ^ 2 + t160) / 0.2e1 + t156 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t160) / 0.2e1 + t154 ^ 2 * Icges(5,3) / 0.2e1 + (m(2) * (t143 ^ 2 + t144 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
