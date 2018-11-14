% Calculate kinetic energy for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:20
% EndTime: 2018-11-14 13:54:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (83->33), mult. (79->55), div. (0->0), fcn. (28->6), ass. (0->26)
t169 = rSges(5,1) + pkin(3);
t159 = qJD(1) + qJD(2);
t168 = pkin(2) * t159;
t167 = pkin(1) * qJD(1);
t160 = qJ(1) + qJ(2);
t162 = cos(qJ(1));
t153 = t162 * t167;
t157 = cos(t160);
t166 = t157 * t168 + t153;
t161 = sin(qJ(1));
t165 = t161 * t167;
t156 = sin(t160);
t164 = -t156 * t168 - t165;
t158 = qJ(3) + t160;
t155 = qJD(3) + t159;
t152 = cos(t158);
t151 = sin(t158);
t150 = t162 * rSges(2,1) - t161 * rSges(2,2);
t149 = t161 * rSges(2,1) + t162 * rSges(2,2);
t147 = t153 + t159 * (t157 * rSges(3,1) - t156 * rSges(3,2));
t146 = -t165 - t159 * (t156 * rSges(3,1) + t157 * rSges(3,2));
t145 = t155 * (t152 * rSges(4,1) - t151 * rSges(4,2)) + t166;
t144 = -t155 * (t151 * rSges(4,1) + t152 * rSges(4,2)) + t164;
t143 = (-t151 * rSges(5,2) + t169 * t152) * t155 + t166;
t142 = (-t152 * rSges(5,2) - t169 * t151) * t155 + t164;
t1 = m(3) * (t146 ^ 2 + t147 ^ 2) / 0.2e1 + t159 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t155 ^ 2 + (m(2) * (t149 ^ 2 + t150 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
