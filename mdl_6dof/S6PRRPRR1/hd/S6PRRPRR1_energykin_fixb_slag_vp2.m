% Calculate kinetic energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:56
% EndTime: 2019-03-08 21:50:57
% DurationCPUTime: 0.45s
% Computational Cost: add. (562->99), mult. (1317->159), div. (0->0), fcn. (989->12), ass. (0->45)
t126 = sin(qJ(5));
t130 = cos(qJ(5));
t121 = sin(pkin(12));
t123 = cos(pkin(12));
t127 = sin(qJ(3));
t131 = cos(qJ(3));
t113 = (t121 * t131 + t123 * t127) * qJD(2);
t128 = sin(qJ(2));
t122 = sin(pkin(6));
t136 = qJD(1) * t122;
t115 = qJD(2) * pkin(8) + t128 * t136;
t124 = cos(pkin(6));
t135 = qJD(1) * t124;
t118 = t131 * t135;
t137 = qJ(4) * qJD(2);
t107 = qJD(3) * pkin(3) + t118 + (-t115 - t137) * t127;
t110 = t131 * t115 + t127 * t135;
t108 = t131 * t137 + t110;
t96 = t123 * t107 - t108 * t121;
t94 = qJD(3) * pkin(4) - pkin(9) * t113 + t96;
t112 = (-t121 * t127 + t123 * t131) * qJD(2);
t97 = t121 * t107 + t123 * t108;
t95 = pkin(9) * t112 + t97;
t90 = t126 * t94 + t130 * t95;
t132 = cos(qJ(2));
t134 = t132 * t136;
t89 = -t126 * t95 + t130 * t94;
t101 = t112 * t130 - t113 * t126;
t111 = -t134 + qJD(4) + (-pkin(3) * t131 - pkin(2)) * qJD(2);
t103 = -pkin(4) * t112 + t111;
t129 = cos(qJ(6));
t125 = sin(qJ(6));
t120 = qJD(3) + qJD(5);
t116 = -qJD(2) * pkin(2) - t134;
t109 = -t115 * t127 + t118;
t102 = t112 * t126 + t113 * t130;
t100 = qJD(6) - t101;
t99 = t102 * t129 + t120 * t125;
t98 = -t102 * t125 + t120 * t129;
t91 = -pkin(5) * t101 - pkin(10) * t102 + t103;
t88 = pkin(10) * t120 + t90;
t87 = -pkin(5) * t120 - t89;
t86 = t125 * t91 + t129 * t88;
t85 = -t125 * t88 + t129 * t91;
t1 = m(4) * (t109 ^ 2 + t110 ^ 2 + t116 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + (t87 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,1) * t99 / 0.2e1) * t99 + (t89 * mrSges(6,1) - t90 * mrSges(6,2) + Ifges(6,3) * t120 / 0.2e1) * t120 + (t111 * mrSges(5,2) - t96 * mrSges(5,3) + Ifges(5,1) * t113 / 0.2e1) * t113 + (-t87 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,2) * t98 / 0.2e1) * t98 + (-t111 * mrSges(5,1) + t97 * mrSges(5,3) + Ifges(5,4) * t113 + Ifges(5,2) * t112 / 0.2e1) * t112 + (t103 * mrSges(6,2) - t89 * mrSges(6,3) + Ifges(6,5) * t120 + Ifges(6,1) * t102 / 0.2e1) * t102 + (m(3) * (t124 ^ 2 + (t128 ^ 2 + t132 ^ 2) * t122 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t103 * mrSges(6,1) + t90 * mrSges(6,3) + Ifges(6,4) * t102 + Ifges(6,6) * t120 + Ifges(6,2) * t101 / 0.2e1) * t101 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,5) * t99 + Ifges(7,6) * t98 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t116 * (-mrSges(4,1) * t131 + mrSges(4,2) * t127) + (Ifges(4,2) * t131 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t131 + Ifges(4,1) * t127 / 0.2e1) * t127) * qJD(2) + (mrSges(3,1) * t132 - mrSges(3,2) * t128) * t136 + (-t109 * t127 + t110 * t131) * mrSges(4,3)) * qJD(2) + (t109 * mrSges(4,1) + t96 * mrSges(5,1) - t110 * mrSges(4,2) - t97 * mrSges(5,2) + Ifges(5,5) * t113 + Ifges(5,6) * t112 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t127 + Ifges(4,6) * t131) * qJD(2)) * qJD(3);
T  = t1;
