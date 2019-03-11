% Calculate kinetic energy for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:53:57
% EndTime: 2019-03-09 06:53:57
% DurationCPUTime: 0.57s
% Computational Cost: add. (594->98), mult. (1263->157), div. (0->0), fcn. (870->10), ass. (0->43)
t125 = m(3) / 0.2e1;
t114 = sin(qJ(5));
t118 = cos(qJ(5));
t115 = sin(qJ(4));
t116 = sin(qJ(3));
t119 = cos(qJ(4));
t120 = cos(qJ(3));
t102 = (t115 * t120 + t116 * t119) * qJD(1);
t110 = qJD(3) + qJD(4);
t111 = sin(pkin(11));
t105 = (pkin(1) * t111 + pkin(7)) * qJD(1);
t109 = t120 * qJD(2);
t124 = pkin(8) * qJD(1);
t97 = qJD(3) * pkin(3) + t109 + (-t105 - t124) * t116;
t100 = t116 * qJD(2) + t120 * t105;
t98 = t120 * t124 + t100;
t86 = -t115 * t98 + t119 * t97;
t83 = pkin(4) * t110 - pkin(9) * t102 + t86;
t101 = (-t115 * t116 + t119 * t120) * qJD(1);
t87 = t115 * t97 + t119 * t98;
t85 = pkin(9) * t101 + t87;
t80 = t114 * t83 + t118 * t85;
t112 = cos(pkin(11));
t123 = -pkin(1) * t112 - pkin(2);
t79 = -t114 * t85 + t118 * t83;
t91 = t101 * t118 - t102 * t114;
t103 = (-pkin(3) * t120 + t123) * qJD(1);
t93 = -pkin(4) * t101 + t103;
t117 = cos(qJ(6));
t113 = sin(qJ(6));
t107 = qJD(5) + t110;
t106 = t123 * qJD(1);
t99 = -t105 * t116 + t109;
t92 = t101 * t114 + t102 * t118;
t90 = qJD(6) - t91;
t89 = t107 * t113 + t117 * t92;
t88 = t107 * t117 - t113 * t92;
t81 = -pkin(5) * t91 - pkin(10) * t92 + t93;
t78 = pkin(10) * t107 + t80;
t77 = -pkin(5) * t107 - t79;
t76 = t113 * t81 + t117 * t78;
t75 = -t113 * t78 + t117 * t81;
t1 = m(7) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t106 ^ 2 + t99 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t125 + (t93 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t92 / 0.2e1) * t92 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,3) * t90 / 0.2e1) * t90 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,3) * t110 / 0.2e1) * t110 + (-t93 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t92 + Ifges(6,2) * t91 / 0.2e1) * t91 + (t77 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,5) * t90 + Ifges(7,1) * t89 / 0.2e1) * t89 + (t103 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,5) * t110 + Ifges(5,1) * t102 / 0.2e1) * t102 + (-t77 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t89 + Ifges(7,6) * t90 + Ifges(7,2) * t88 / 0.2e1) * t88 + (t79 * mrSges(6,1) - t80 * mrSges(6,2) + Ifges(6,5) * t92 + Ifges(6,6) * t91 + Ifges(6,3) * t107 / 0.2e1) * t107 + (-t103 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,4) * t102 + Ifges(5,6) * t110 + Ifges(5,2) * t101 / 0.2e1) * t101 + (t99 * mrSges(4,1) - t100 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t106 * (-mrSges(4,1) * t120 + mrSges(4,2) * t116) + (t100 * t120 - t99 * t116) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t116 + Ifges(4,6) * t120) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t112 * mrSges(3,1) - t111 * mrSges(3,2) + (t111 ^ 2 + t112 ^ 2) * t125 * pkin(1)) * pkin(1) + Ifges(4,2) * t120 ^ 2 / 0.2e1 + (Ifges(4,4) * t120 + Ifges(4,1) * t116 / 0.2e1) * t116) * qJD(1)) * qJD(1);
T  = t1;
