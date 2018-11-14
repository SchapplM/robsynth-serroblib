% Calculate kinetic energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:03
% EndTime: 2018-11-14 13:38:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->18), mult. (47->34), div. (0->0), fcn. (26->4), ass. (0->12)
t147 = sin(pkin(5));
t148 = cos(pkin(5));
t142 = qJD(2) * t147 + qJD(3) * t148;
t143 = -qJD(2) * t148 + qJD(3) * t147;
t152 = qJD(1) ^ 2;
t150 = cos(qJ(4));
t149 = sin(qJ(4));
t141 = t147 * t150 + t148 * t149;
t140 = -t147 * t149 + t148 * t150;
t139 = qJD(4) * (t141 * rSges(5,1) + t140 * rSges(5,2)) + t143;
t138 = -qJD(4) * (-t140 * rSges(5,1) + t141 * rSges(5,2)) + t142;
t1 = m(2) * t152 / 0.2e1 + m(3) * (t152 + (t147 ^ 2 + t148 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t152) / 0.2e1 + m(5) * (t138 ^ 2 + t139 ^ 2 + t152) / 0.2e1 + qJD(4) ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
