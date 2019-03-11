% Calculate kinetic energy for
% S4PPRR2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:49
% EndTime: 2019-03-08 18:16:49
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->21), mult. (40->35), div. (0->0), fcn. (10->4), ass. (0->16)
t125 = pkin(3) * qJD(3);
t120 = pkin(6) + qJ(3);
t124 = qJD(1) ^ 2;
t123 = qJD(2) ^ 2;
t122 = qJD(3) ^ 2;
t121 = qJD(3) + qJD(4);
t119 = qJ(4) + t120;
t118 = cos(t120);
t117 = sin(t120);
t116 = cos(t119);
t115 = sin(t119);
t114 = t117 * rSges(4,1) + t118 * rSges(4,2);
t113 = qJD(1) + qJD(3) * (t118 * rSges(4,1) - t117 * rSges(4,2));
t112 = -t117 * t125 - t121 * (t115 * rSges(5,1) + t116 * rSges(5,2));
t111 = qJD(1) + t118 * t125 + t121 * (t116 * rSges(5,1) - t115 * rSges(5,2));
t1 = m(2) * t124 / 0.2e1 + m(3) * (t123 + t124) / 0.2e1 + m(4) * (t122 * t114 ^ 2 + t113 ^ 2 + t123) / 0.2e1 + t122 * Icges(4,3) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t123) / 0.2e1 + t121 ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
