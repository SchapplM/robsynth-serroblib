% Calculate kinetic energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:57
% EndTime: 2019-03-09 06:56:57
% DurationCPUTime: 0.58s
% Computational Cost: add. (592->98), mult. (1185->158), div. (0->0), fcn. (798->10), ass. (0->43)
t129 = m(3) / 0.2e1;
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t114 = qJD(3) + qJD(4);
t115 = sin(pkin(11));
t109 = (pkin(1) * t115 + pkin(7)) * qJD(1);
t120 = sin(qJ(3));
t124 = cos(qJ(3));
t102 = t120 * qJD(2) + t124 * t109;
t128 = qJD(1) * t124;
t100 = pkin(8) * t128 + t102;
t119 = sin(qJ(4));
t123 = cos(qJ(4));
t113 = t124 * qJD(2);
t97 = qJD(3) * pkin(3) + t113 + (-pkin(8) * qJD(1) - t109) * t120;
t90 = t123 * t100 + t119 * t97;
t86 = pkin(9) * t114 + t90;
t105 = -t119 * t120 * qJD(1) + t123 * t128;
t106 = (t119 * t124 + t120 * t123) * qJD(1);
t116 = cos(pkin(11));
t127 = -pkin(1) * t116 - pkin(2);
t107 = (-pkin(3) * t124 + t127) * qJD(1);
t93 = -pkin(4) * t105 - pkin(9) * t106 + t107;
t82 = t118 * t93 + t122 * t86;
t81 = -t118 * t86 + t122 * t93;
t89 = -t119 * t100 + t123 * t97;
t85 = -pkin(4) * t114 - t89;
t104 = qJD(5) - t105;
t121 = cos(qJ(6));
t117 = sin(qJ(6));
t110 = t127 * qJD(1);
t103 = qJD(6) + t104;
t101 = -t109 * t120 + t113;
t99 = t106 * t122 + t114 * t118;
t98 = -t106 * t118 + t114 * t122;
t88 = t117 * t98 + t121 * t99;
t87 = -t117 * t99 + t121 * t98;
t83 = -pkin(5) * t98 + t85;
t80 = pkin(10) * t98 + t82;
t79 = pkin(5) * t104 - pkin(10) * t99 + t81;
t78 = t117 * t79 + t121 * t80;
t77 = -t117 * t80 + t121 * t79;
t1 = m(7) * (t77 ^ 2 + t78 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t129 + (t85 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,1) * t99 / 0.2e1) * t99 + (t83 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t88 / 0.2e1) * t88 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (-t85 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,2) * t98 / 0.2e1) * t98 + (-t83 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t107 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t106 / 0.2e1) * t106 + (-t107 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,6) * t114 + Ifges(5,2) * t105 / 0.2e1) * t105 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,5) * t99 + Ifges(6,6) * t98 + Ifges(6,3) * t104 / 0.2e1) * t104 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t103 / 0.2e1) * t103 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t110 * (-mrSges(4,1) * t124 + mrSges(4,2) * t120) + (-t101 * t120 + t102 * t124) * mrSges(4,3) + (Ifges(4,5) * t120 + Ifges(4,6) * t124) * qJD(3) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t116 * mrSges(3,1) - t115 * mrSges(3,2) + (t115 ^ 2 + t116 ^ 2) * t129 * pkin(1)) * pkin(1) + Ifges(4,2) * t124 ^ 2 / 0.2e1 + (Ifges(4,4) * t124 + Ifges(4,1) * t120 / 0.2e1) * t120) * qJD(1)) * qJD(1);
T  = t1;
