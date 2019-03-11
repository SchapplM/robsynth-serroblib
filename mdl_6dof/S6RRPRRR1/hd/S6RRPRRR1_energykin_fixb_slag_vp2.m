% Calculate kinetic energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:23
% EndTime: 2019-03-09 13:12:24
% DurationCPUTime: 0.71s
% Computational Cost: add. (1016->110), mult. (2495->171), div. (0->0), fcn. (1926->10), ass. (0->44)
t129 = qJD(1) * (pkin(7) + qJ(3));
t127 = pkin(7) * mrSges(3,3);
t117 = sin(qJ(5));
t121 = cos(qJ(5));
t114 = sin(pkin(11));
t115 = cos(pkin(11));
t119 = sin(qJ(2));
t123 = cos(qJ(2));
t106 = (-t114 * t119 + t115 * t123) * qJD(1);
t107 = (t114 * t123 + t115 * t119) * qJD(1);
t118 = sin(qJ(4));
t122 = cos(qJ(4));
t100 = t106 * t118 + t107 * t122;
t113 = qJD(2) + qJD(4);
t109 = qJD(2) * pkin(2) - t119 * t129;
t110 = t123 * t129;
t101 = t115 * t109 - t110 * t114;
t97 = qJD(2) * pkin(3) - pkin(8) * t107 + t101;
t102 = t114 * t109 + t115 * t110;
t98 = pkin(8) * t106 + t102;
t86 = -t118 * t98 + t122 * t97;
t83 = pkin(4) * t113 - pkin(9) * t100 + t86;
t87 = t118 * t97 + t122 * t98;
t99 = t106 * t122 - t107 * t118;
t85 = pkin(9) * t99 + t87;
t80 = t117 * t83 + t121 * t85;
t79 = -t117 * t85 + t121 * t83;
t91 = -t100 * t117 + t121 * t99;
t111 = qJD(3) + (-pkin(2) * t123 - pkin(1)) * qJD(1);
t103 = -pkin(3) * t106 + t111;
t93 = -pkin(4) * t99 + t103;
t120 = cos(qJ(6));
t116 = sin(qJ(6));
t112 = qJD(5) + t113;
t92 = t100 * t121 + t117 * t99;
t90 = qJD(6) - t91;
t89 = t112 * t116 + t120 * t92;
t88 = t112 * t120 - t116 * t92;
t81 = -pkin(5) * t91 - pkin(10) * t92 + t93;
t78 = pkin(10) * t112 + t80;
t77 = -pkin(5) * t112 - t79;
t76 = t116 * t81 + t120 * t78;
t75 = -t116 * t78 + t120 * t81;
t1 = m(7) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t93 ^ 2) / 0.2e1 + (-t103 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,2) * t99 / 0.2e1) * t99 + (t93 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t92 / 0.2e1) * t92 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,3) * t90 / 0.2e1) * t90 + (t111 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,1) * t107 / 0.2e1) * t107 + (-t93 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t92 + Ifges(6,2) * t91 / 0.2e1) * t91 + (t77 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,5) * t90 + Ifges(7,1) * t89 / 0.2e1) * t89 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,6) * t99 + Ifges(5,3) * t113 / 0.2e1) * t113 + (-t111 * mrSges(4,1) + t102 * mrSges(4,3) + Ifges(4,4) * t107 + Ifges(4,2) * t106 / 0.2e1) * t106 + (-t77 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t89 + Ifges(7,6) * t90 + Ifges(7,2) * t88 / 0.2e1) * t88 + (t79 * mrSges(6,1) - t80 * mrSges(6,2) + Ifges(6,5) * t92 + Ifges(6,6) * t91 + Ifges(6,3) * t112 / 0.2e1) * t112 + (t103 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,4) * t99 + Ifges(5,5) * t113 + Ifges(5,1) * t100 / 0.2e1) * t100 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,5) * t107 + Ifges(4,6) * t106 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t119 + Ifges(3,6) * t123 + (-mrSges(3,1) * t119 - mrSges(3,2) * t123) * pkin(7)) * qJD(1)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t119 ^ 2 + t123 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t127 + Ifges(3,2) / 0.2e1) * t123) * t123 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t123 + (t127 + Ifges(3,1) / 0.2e1) * t119) * t119) * qJD(1) ^ 2;
T  = t1;
