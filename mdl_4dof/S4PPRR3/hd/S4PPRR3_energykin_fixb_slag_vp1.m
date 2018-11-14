% Calculate kinetic energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:15
% EndTime: 2018-11-14 14:08:15
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->21), mult. (40->35), div. (0->0), fcn. (10->4), ass. (0->16)
t124 = pkin(3) * qJD(3);
t119 = pkin(6) + qJ(3);
t123 = qJD(1) ^ 2;
t122 = qJD(2) ^ 2;
t121 = qJD(3) ^ 2;
t120 = -qJD(3) - qJD(4);
t118 = qJ(4) + t119;
t117 = cos(t119);
t116 = sin(t119);
t115 = cos(t118);
t114 = sin(t118);
t113 = t116 * rSges(4,1) + t117 * rSges(4,2);
t112 = qJD(1) + qJD(3) * (t117 * rSges(4,1) - t116 * rSges(4,2));
t111 = -t116 * t124 + t120 * (t114 * rSges(5,1) + t115 * rSges(5,2));
t110 = qJD(1) + t117 * t124 - t120 * (t115 * rSges(5,1) - t114 * rSges(5,2));
t1 = m(2) * t123 / 0.2e1 + m(3) * (t122 + t123) / 0.2e1 + m(4) * (t121 * t113 ^ 2 + t112 ^ 2 + t122) / 0.2e1 + t121 * Icges(4,3) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t122) / 0.2e1 + t120 ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
