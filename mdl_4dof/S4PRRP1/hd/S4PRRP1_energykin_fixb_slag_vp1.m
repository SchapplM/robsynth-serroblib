% Calculate kinetic energy for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:22
% EndTime: 2018-11-14 13:43:22
% DurationCPUTime: 0.04s
% Computational Cost: add. (72->29), mult. (63->44), div. (0->0), fcn. (22->4), ass. (0->21)
t165 = rSges(5,1) + pkin(3);
t164 = pkin(2) * qJD(2);
t163 = rSges(5,3) + qJ(4);
t158 = pkin(6) + qJ(2);
t154 = sin(t158);
t162 = t154 * t164;
t161 = qJD(1) ^ 2;
t160 = qJD(2) ^ 2;
t159 = qJD(2) + qJD(3);
t156 = qJ(3) + t158;
t155 = cos(t158);
t153 = cos(t156);
t152 = sin(t156);
t151 = t155 * t164;
t150 = t155 * rSges(3,1) - t154 * rSges(3,2);
t149 = t154 * rSges(3,1) + t155 * rSges(3,2);
t148 = t151 + t159 * (t153 * rSges(4,1) - t152 * rSges(4,2));
t147 = -t162 - t159 * (t152 * rSges(4,1) + t153 * rSges(4,2));
t146 = -qJD(4) * t153 + t151 + (t163 * t152 + t165 * t153) * t159;
t145 = -t162 + qJD(4) * t152 + (-t165 * t152 + t163 * t153) * t159;
t1 = m(2) * t161 / 0.2e1 + m(3) * (t161 + (t149 ^ 2 + t150 ^ 2) * t160) / 0.2e1 + t160 * Icges(3,3) / 0.2e1 + m(4) * (t147 ^ 2 + t148 ^ 2 + t161) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t161) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t159 ^ 2;
T  = t1;
