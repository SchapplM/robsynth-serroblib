% Calculate kinetic energy for
% S6RRRRRR1
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:48
% EndTime: 2019-03-10 03:27:48
% DurationCPUTime: 0.72s
% Computational Cost: add. (1076->110), mult. (2495->172), div. (0->0), fcn. (1926->10), ass. (0->45)
t130 = qJD(1) * (-pkin(8) - pkin(7));
t128 = pkin(7) * mrSges(3,3);
t117 = sin(qJ(5));
t122 = cos(qJ(5));
t119 = sin(qJ(3));
t120 = sin(qJ(2));
t124 = cos(qJ(3));
t125 = cos(qJ(2));
t107 = (-t119 * t120 + t124 * t125) * qJD(1);
t108 = (t119 * t125 + t120 * t124) * qJD(1);
t118 = sin(qJ(4));
t123 = cos(qJ(4));
t101 = t107 * t118 + t108 * t123;
t115 = qJD(2) + qJD(3);
t114 = qJD(4) + t115;
t110 = qJD(2) * pkin(2) + t120 * t130;
t111 = t125 * t130;
t102 = t124 * t110 + t111 * t119;
t97 = pkin(3) * t115 - pkin(9) * t108 + t102;
t103 = t119 * t110 - t124 * t111;
t99 = pkin(9) * t107 + t103;
t87 = -t118 * t99 + t123 * t97;
t84 = pkin(4) * t114 - pkin(10) * t101 + t87;
t100 = t107 * t123 - t108 * t118;
t88 = t118 * t97 + t123 * t99;
t86 = pkin(10) * t100 + t88;
t81 = t117 * t84 + t122 * t86;
t112 = (-pkin(2) * t125 - pkin(1)) * qJD(1);
t80 = -t117 * t86 + t122 * t84;
t92 = t100 * t122 - t101 * t117;
t104 = -pkin(3) * t107 + t112;
t94 = -pkin(4) * t100 + t104;
t121 = cos(qJ(6));
t116 = sin(qJ(6));
t113 = qJD(5) + t114;
t93 = t100 * t117 + t101 * t122;
t91 = qJD(6) - t92;
t90 = t113 * t116 + t121 * t93;
t89 = t113 * t121 - t116 * t93;
t82 = -pkin(5) * t92 - pkin(11) * t93 + t94;
t79 = pkin(11) * t113 + t81;
t78 = -pkin(5) * t113 - t80;
t77 = t116 * t82 + t121 * t79;
t76 = -t116 * t79 + t121 * t82;
t1 = m(4) * (t102 ^ 2 + t103 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t94 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,3) * t91 / 0.2e1) * t91 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,3) * t115 / 0.2e1) * t115 + (t87 * mrSges(5,1) - t88 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (-t94 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (t78 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,5) * t91 + Ifges(7,1) * t90 / 0.2e1) * t90 + (t112 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,5) * t115 + Ifges(4,1) * t108 / 0.2e1) * t108 + (t104 * mrSges(5,2) - t87 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t101 / 0.2e1) * t101 + (-t78 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t90 + Ifges(7,6) * t91 + Ifges(7,2) * t89 / 0.2e1) * t89 + (t80 * mrSges(6,1) - t81 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t113 / 0.2e1) * t113 + (-t112 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t108 + Ifges(4,6) * t115 + Ifges(4,2) * t107 / 0.2e1) * t107 + (-t104 * mrSges(5,1) + t88 * mrSges(5,3) + Ifges(5,4) * t101 + Ifges(5,6) * t114 + Ifges(5,2) * t100 / 0.2e1) * t100 + ((m(3) * (pkin(1) ^ 2 + (t120 ^ 2 + t125 ^ 2) * pkin(7) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t128) * t125) * t125 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t125 + (Ifges(3,1) / 0.2e1 + t128) * t120) * t120) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t125 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t120) * qJD(2)) * qJD(1);
T  = t1;
