% Calculate kinetic energy for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:27
% EndTime: 2019-03-09 00:37:27
% DurationCPUTime: 0.48s
% Computational Cost: add. (594->99), mult. (1317->160), div. (0->0), fcn. (989->12), ass. (0->46)
t126 = sin(qJ(5));
t131 = cos(qJ(5));
t127 = sin(qJ(4));
t128 = sin(qJ(3));
t132 = cos(qJ(4));
t133 = cos(qJ(3));
t114 = (t127 * t133 + t128 * t132) * qJD(2);
t122 = qJD(3) + qJD(4);
t129 = sin(qJ(2));
t123 = sin(pkin(6));
t139 = qJD(1) * t123;
t116 = qJD(2) * pkin(8) + t129 * t139;
t124 = cos(pkin(6));
t138 = qJD(1) * t124;
t119 = t133 * t138;
t108 = qJD(3) * pkin(3) + t119 + (-pkin(9) * qJD(2) - t116) * t128;
t111 = t133 * t116 + t128 * t138;
t137 = qJD(2) * t133;
t109 = pkin(9) * t137 + t111;
t97 = t132 * t108 - t109 * t127;
t94 = pkin(4) * t122 - pkin(10) * t114 + t97;
t113 = (-t127 * t128 + t132 * t133) * qJD(2);
t98 = t127 * t108 + t132 * t109;
t96 = pkin(10) * t113 + t98;
t91 = t126 * t94 + t131 * t96;
t134 = cos(qJ(2));
t136 = t134 * t139;
t90 = -t126 * t96 + t131 * t94;
t102 = t113 * t131 - t114 * t126;
t112 = -t136 + (-pkin(3) * t133 - pkin(2)) * qJD(2);
t106 = -pkin(4) * t113 + t112;
t130 = cos(qJ(6));
t125 = sin(qJ(6));
t120 = qJD(5) + t122;
t117 = -qJD(2) * pkin(2) - t136;
t110 = -t116 * t128 + t119;
t103 = t113 * t126 + t114 * t131;
t101 = qJD(6) - t102;
t100 = t103 * t130 + t120 * t125;
t99 = -t103 * t125 + t120 * t130;
t92 = -pkin(5) * t102 - pkin(11) * t103 + t106;
t89 = pkin(11) * t120 + t91;
t88 = -pkin(5) * t120 - t90;
t87 = t125 * t92 + t130 * t89;
t86 = -t125 * t89 + t130 * t92;
t1 = m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t117 ^ 2) / 0.2e1 + (-t88 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t97 * mrSges(5,1) - t98 * mrSges(5,2) + Ifges(5,3) * t122 / 0.2e1) * t122 + (t90 * mrSges(6,1) - t91 * mrSges(6,2) + Ifges(6,3) * t120 / 0.2e1) * t120 + (t110 * mrSges(4,1) - t111 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t112 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,5) * t122 + Ifges(5,1) * t114 / 0.2e1) * t114 + (t106 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,5) * t120 + Ifges(6,1) * t103 / 0.2e1) * t103 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t101 / 0.2e1) * t101 + (m(3) * (t124 ^ 2 + (t129 ^ 2 + t134 ^ 2) * t123 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t112 * mrSges(5,1) + t98 * mrSges(5,3) + Ifges(5,4) * t114 + Ifges(5,6) * t122 + Ifges(5,2) * t113 / 0.2e1) * t113 + (-t106 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t103 + Ifges(6,6) * t120 + Ifges(6,2) * t102 / 0.2e1) * t102 + (t88 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t101 + Ifges(7,1) * t100 / 0.2e1) * t100 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t134 - mrSges(3,2) * t129) * t139 + (-t117 * mrSges(4,1) + t111 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t137 / 0.2e1) * t133 + (t117 * mrSges(4,2) - t110 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t133 + Ifges(4,1) * t128 / 0.2e1) * qJD(2)) * t128) * qJD(2);
T  = t1;
