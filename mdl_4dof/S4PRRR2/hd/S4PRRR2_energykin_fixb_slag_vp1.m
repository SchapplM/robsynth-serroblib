% Calculate kinetic energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:21
% EndTime: 2019-07-18 13:27:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (49->28), mult. (59->47), div. (0->0), fcn. (18->6), ass. (0->23)
t111 = -qJD(2) - qJD(3);
t120 = pkin(2) * t111;
t119 = pkin(1) * qJD(2);
t112 = qJ(2) + qJ(3);
t113 = sin(qJ(2));
t118 = t113 * t119;
t114 = cos(qJ(2));
t117 = t114 * t119;
t116 = qJD(1) ^ 2;
t115 = qJD(2) ^ 2;
t110 = qJ(4) + t112;
t109 = cos(t112);
t108 = sin(t112);
t107 = -qJD(4) + t111;
t106 = cos(t110);
t105 = sin(t110);
t104 = t114 * rSges(3,1) - t113 * rSges(3,2);
t103 = -t113 * rSges(3,1) - t114 * rSges(3,2);
t102 = -t117 + t111 * (t109 * rSges(4,1) - t108 * rSges(4,2));
t101 = -t118 - t111 * (-t108 * rSges(4,1) - t109 * rSges(4,2));
t100 = -t117 + t109 * t120 + t107 * (t106 * rSges(5,1) - t105 * rSges(5,2));
t99 = -t118 + t108 * t120 - t107 * (-t105 * rSges(5,1) - t106 * rSges(5,2));
t1 = m(2) * t116 / 0.2e1 + m(3) * (t116 + (t103 ^ 2 + t104 ^ 2) * t115) / 0.2e1 + t115 * Icges(3,3) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t116) / 0.2e1 + t111 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t100 ^ 2 + t99 ^ 2 + t116) / 0.2e1 + t107 ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
