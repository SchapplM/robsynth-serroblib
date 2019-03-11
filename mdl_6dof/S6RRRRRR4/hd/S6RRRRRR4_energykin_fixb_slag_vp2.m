% Calculate kinetic energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:38
% EndTime: 2019-03-10 03:46:38
% DurationCPUTime: 0.69s
% Computational Cost: add. (1178->111), mult. (2493->173), div. (0->0), fcn. (1874->10), ass. (0->46)
t133 = pkin(7) * mrSges(3,3);
t120 = sin(qJ(5));
t125 = cos(qJ(5));
t122 = sin(qJ(3));
t127 = cos(qJ(3));
t123 = sin(qJ(2));
t132 = qJD(1) * t123;
t109 = qJD(2) * t127 - t122 * t132;
t110 = qJD(2) * t122 + t127 * t132;
t121 = sin(qJ(4));
t126 = cos(qJ(4));
t102 = t109 * t121 + t110 * t126;
t128 = cos(qJ(2));
t131 = qJD(1) * t128;
t117 = qJD(3) - t131;
t116 = qJD(4) + t117;
t108 = (-pkin(2) * t128 - pkin(8) * t123 - pkin(1)) * qJD(1);
t115 = pkin(7) * t131 + qJD(2) * pkin(8);
t104 = t108 * t122 + t115 * t127;
t100 = pkin(9) * t109 + t104;
t103 = t108 * t127 - t115 * t122;
t98 = pkin(3) * t117 - pkin(9) * t110 + t103;
t91 = -t100 * t121 + t126 * t98;
t87 = pkin(4) * t116 - pkin(10) * t102 + t91;
t101 = t109 * t126 - t110 * t121;
t92 = t100 * t126 + t121 * t98;
t89 = pkin(10) * t101 + t92;
t82 = t120 * t87 + t125 * t89;
t81 = -t120 * t89 + t125 * t87;
t114 = -qJD(2) * pkin(2) + pkin(7) * t132;
t105 = -pkin(3) * t109 + t114;
t113 = qJD(5) + t116;
t95 = -pkin(4) * t101 + t105;
t124 = cos(qJ(6));
t119 = sin(qJ(6));
t112 = qJD(6) + t113;
t94 = t101 * t120 + t102 * t125;
t93 = t101 * t125 - t102 * t120;
t90 = -pkin(5) * t93 + t95;
t86 = t119 * t93 + t124 * t94;
t85 = -t119 * t94 + t124 * t93;
t80 = pkin(11) * t93 + t82;
t79 = pkin(5) * t113 - pkin(11) * t94 + t81;
t78 = t119 * t79 + t124 * t80;
t77 = -t119 * t80 + t124 * t79;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(5) * (t105 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t114 ^ 2) / 0.2e1 + (t95 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,1) * t94 / 0.2e1) * t94 + (t90 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t86 / 0.2e1) * t86 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,3) * t117 / 0.2e1) * t117 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t116 / 0.2e1) * t116 + (-t95 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t94 + Ifges(6,2) * t93 / 0.2e1) * t93 + (-t90 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t114 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,5) * t117 + Ifges(4,1) * t110 / 0.2e1) * t110 + (t105 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t116 + Ifges(5,1) * t102 / 0.2e1) * t102 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,5) * t94 + Ifges(6,6) * t93 + Ifges(6,3) * t113 / 0.2e1) * t113 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t86 + Ifges(7,6) * t85 + Ifges(7,3) * t112 / 0.2e1) * t112 + (-t114 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t110 + Ifges(4,6) * t117 + Ifges(4,2) * t109 / 0.2e1) * t109 + (-t105 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t102 + Ifges(5,6) * t116 + Ifges(5,2) * t101 / 0.2e1) * t101 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t123 ^ 2 + t128 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t133 + Ifges(3,2) / 0.2e1) * t128) * t128 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t128 + (t133 + Ifges(3,1) / 0.2e1) * t123) * t123) * qJD(1) + ((-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t128 + (-mrSges(3,1) * pkin(7) + Ifges(3,5)) * t123) * qJD(2)) * qJD(1);
T  = t1;
