% Calculate kinetic energy for
% S4PRPR1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->34), mult. (75->53), div. (0->0), fcn. (38->4), ass. (0->19)
t160 = qJD(1) ^ 2;
t158 = cos(qJ(4));
t157 = sin(qJ(4));
t156 = qJD(2) - qJD(4);
t155 = pkin(6) + qJ(2);
t154 = cos(t155);
t153 = sin(t155);
t152 = qJD(3) * t153;
t151 = t154 * rSges(3,1) - t153 * rSges(3,2);
t150 = t153 * rSges(3,1) + t154 * rSges(3,2);
t149 = t153 * pkin(2) - t154 * qJ(3);
t148 = qJD(2) * (t154 * pkin(2) + t153 * qJ(3));
t147 = t153 * t158 - t154 * t157;
t146 = -t153 * t157 - t154 * t158;
t145 = t148 - qJD(3) * t154 + qJD(2) * (t154 * rSges(4,1) + t153 * rSges(4,3));
t144 = t152 + (-t153 * rSges(4,1) + t154 * rSges(4,3) - t149) * qJD(2);
t143 = t148 + t156 * (-t146 * rSges(5,1) + t147 * rSges(5,2)) + (qJD(2) * pkin(3) - qJD(3)) * t154;
t142 = t152 - t156 * (t147 * rSges(5,1) + t146 * rSges(5,2)) + (-t153 * pkin(3) - t149) * qJD(2);
t1 = m(4) * (t144 ^ 2 + t145 ^ 2 + t160) / 0.2e1 + m(5) * (t142 ^ 2 + t143 ^ 2 + t160) / 0.2e1 + t156 ^ 2 * Icges(5,3) / 0.2e1 + (m(3) * (t150 ^ 2 + t151 ^ 2) / 0.2e1 + Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * qJD(2) ^ 2 + (m(2) + m(3)) * t160 / 0.2e1;
T  = t1;
