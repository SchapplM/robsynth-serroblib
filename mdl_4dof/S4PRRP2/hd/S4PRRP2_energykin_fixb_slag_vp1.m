% Calculate kinetic energy for
% S4PRRP2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:59
% EndTime: 2019-03-08 18:23:59
% DurationCPUTime: 0.03s
% Computational Cost: add. (43->25), mult. (56->44), div. (0->0), fcn. (18->4), ass. (0->18)
t124 = rSges(5,1) + pkin(3);
t123 = pkin(2) * qJD(2);
t119 = cos(qJ(2));
t122 = t119 * t123 + qJD(1);
t118 = sin(qJ(2));
t121 = t118 * t123;
t120 = qJD(2) ^ 2;
t117 = qJ(2) + qJ(3);
t116 = qJD(2) + qJD(3);
t114 = cos(t117);
t113 = sin(t117);
t111 = t118 * rSges(3,1) + t119 * rSges(3,2);
t110 = qJD(1) + qJD(2) * (t119 * rSges(3,1) - t118 * rSges(3,2));
t109 = -t121 - t116 * (t113 * rSges(4,1) + t114 * rSges(4,2));
t108 = t116 * (t114 * rSges(4,1) - t113 * rSges(4,2)) + t122;
t107 = -t121 + (-t114 * rSges(5,2) - t124 * t113) * t116;
t106 = (-t113 * rSges(5,2) + t124 * t114) * t116 + t122;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t120 * t111 ^ 2 + t110 ^ 2) / 0.2e1 + t120 * Icges(3,3) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t116 ^ 2;
T  = t1;
