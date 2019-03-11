% Calculate kinetic energy for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:08:59
% EndTime: 2019-03-09 05:09:00
% DurationCPUTime: 0.69s
% Computational Cost: add. (909->104), mult. (2211->162), div. (0->0), fcn. (1724->10), ass. (0->45)
t122 = cos(pkin(10));
t135 = t122 ^ 2;
t134 = qJD(1) * (pkin(7) + qJ(2));
t133 = m(3) / 0.2e1;
t132 = cos(qJ(4));
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t118 = qJD(3) + qJD(4);
t124 = sin(qJ(4));
t120 = sin(pkin(10));
t112 = t120 * t134;
t113 = t122 * t134;
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t103 = -t127 * t112 - t113 * t125;
t111 = (t120 * t127 + t122 * t125) * qJD(1);
t96 = qJD(3) * pkin(3) - pkin(8) * t111 + t103;
t104 = -t125 * t112 + t127 * t113;
t110 = (-t120 * t125 + t122 * t127) * qJD(1);
t99 = pkin(8) * t110 + t104;
t89 = t124 * t96 + t132 * t99;
t85 = qJ(5) * t118 + t89;
t101 = -t132 * t110 + t111 * t124;
t102 = t124 * t110 + t132 * t111;
t114 = qJD(2) + (-pkin(2) * t122 - pkin(1)) * qJD(1);
t105 = -pkin(3) * t110 + t114;
t90 = pkin(4) * t101 - qJ(5) * t102 + t105;
t81 = t119 * t90 + t121 * t85;
t80 = -t119 * t85 + t121 * t90;
t88 = -t124 * t99 + t132 * t96;
t84 = -t118 * pkin(4) + qJD(5) - t88;
t126 = cos(qJ(6));
t123 = sin(qJ(6));
t115 = -qJD(1) * pkin(1) + qJD(2);
t100 = qJD(6) + t101;
t98 = t102 * t121 + t118 * t119;
t97 = -t102 * t119 + t118 * t121;
t92 = t123 * t97 + t126 * t98;
t91 = -t123 * t98 + t126 * t97;
t82 = -t97 * pkin(5) + t84;
t79 = pkin(9) * t97 + t81;
t78 = pkin(5) * t101 - pkin(9) * t98 + t80;
t77 = t123 * t78 + t126 * t79;
t76 = -t123 * t79 + t126 * t78;
t1 = m(7) * (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + t115 ^ 2 * t133 + m(4) * (t103 ^ 2 + t104 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + (t84 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t98 / 0.2e1) * t98 + (t82 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + Ifges(5,3) * t118 / 0.2e1) * t118 + (t114 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,1) * t111 / 0.2e1) * t111 + (-t84 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t98 + Ifges(6,2) * t97 / 0.2e1) * t97 + (-t82 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (-t114 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t111 + Ifges(4,2) * t110 / 0.2e1) * t110 + (t105 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,5) * t118 + Ifges(5,1) * t102 / 0.2e1) * t102 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,5) * t111 + Ifges(4,6) * t110 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t105 * mrSges(5,1) + t80 * mrSges(6,1) - t81 * mrSges(6,2) - t89 * mrSges(5,3) - Ifges(5,4) * t102 + Ifges(6,5) * t98 - Ifges(5,6) * t118 + Ifges(6,6) * t97 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t101) * t101 + (t115 * (-mrSges(3,1) * t122 + mrSges(3,2) * t120) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t133 + mrSges(3,3)) * (t120 ^ 2 + t135) * qJ(2) + Ifges(3,2) * t135 / 0.2e1 + (Ifges(3,4) * t122 + Ifges(3,1) * t120 / 0.2e1) * t120) * qJD(1)) * qJD(1);
T  = t1;
