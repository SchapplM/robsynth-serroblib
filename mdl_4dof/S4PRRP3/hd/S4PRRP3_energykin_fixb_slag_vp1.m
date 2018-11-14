% Calculate kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:13
% EndTime: 2018-11-14 14:12:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (43->25), mult. (56->44), div. (0->0), fcn. (18->4), ass. (0->18)
t125 = rSges(5,1) + pkin(3);
t124 = pkin(2) * qJD(2);
t120 = cos(qJ(2));
t123 = t120 * t124 + qJD(1);
t119 = sin(qJ(2));
t122 = t119 * t124;
t121 = qJD(2) ^ 2;
t118 = qJ(2) + qJ(3);
t117 = -qJD(2) - qJD(3);
t115 = cos(t118);
t114 = sin(t118);
t112 = t119 * rSges(3,1) + t120 * rSges(3,2);
t111 = qJD(1) + qJD(2) * (t120 * rSges(3,1) - t119 * rSges(3,2));
t110 = -t122 + t117 * (t114 * rSges(4,1) + t115 * rSges(4,2));
t109 = -t117 * (t115 * rSges(4,1) - t114 * rSges(4,2)) + t123;
t108 = -t122 + (t115 * rSges(5,2) + t125 * t114) * t117;
t107 = (t114 * rSges(5,2) - t125 * t115) * t117 + t123;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t121 * t112 ^ 2 + t111 ^ 2) / 0.2e1 + t121 * Icges(3,3) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t117 ^ 2;
T  = t1;
