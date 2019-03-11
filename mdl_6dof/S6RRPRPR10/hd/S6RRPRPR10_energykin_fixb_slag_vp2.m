% Calculate kinetic energy for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:18
% EndTime: 2019-03-09 11:04:19
% DurationCPUTime: 0.56s
% Computational Cost: add. (914->111), mult. (2134->163), div. (0->0), fcn. (1666->10), ass. (0->48)
t141 = pkin(4) + pkin(10);
t140 = cos(qJ(4));
t128 = sin(qJ(4));
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t124 = sin(pkin(6));
t138 = t124 * qJD(1);
t135 = t131 * t138;
t139 = cos(pkin(6)) * qJD(1);
t137 = pkin(1) * t139;
t117 = pkin(8) * t135 + t129 * t137;
t122 = qJD(2) + t139;
t111 = qJ(3) * t122 + t117;
t113 = (-pkin(2) * t131 - qJ(3) * t129 - pkin(1)) * t138;
t123 = sin(pkin(11));
t125 = cos(pkin(11));
t101 = -t111 * t123 + t125 * t113;
t136 = t129 * t138;
t115 = t122 * t123 + t125 * t136;
t95 = -pkin(3) * t135 - pkin(9) * t115 + t101;
t102 = t125 * t111 + t123 * t113;
t114 = t122 * t125 - t123 * t136;
t98 = pkin(9) * t114 + t102;
t92 = t128 * t95 + t140 * t98;
t118 = qJD(4) - t135;
t90 = -qJ(5) * t118 - t92;
t91 = -t128 * t98 + t140 * t95;
t116 = -pkin(8) * t136 + t131 * t137;
t134 = qJD(5) - t91;
t106 = t128 * t114 + t140 * t115;
t108 = -pkin(2) * t122 + qJD(3) - t116;
t104 = -pkin(3) * t114 + t108;
t133 = -qJ(5) * t106 + t104;
t132 = qJD(1) ^ 2;
t130 = cos(qJ(6));
t127 = sin(qJ(6));
t105 = -t140 * t114 + t115 * t128;
t103 = qJD(6) + t106;
t100 = t105 * t127 + t118 * t130;
t99 = t105 * t130 - t118 * t127;
t93 = pkin(4) * t105 + t133;
t89 = -t118 * pkin(4) + t134;
t88 = t141 * t105 + t133;
t87 = -pkin(5) * t105 - t90;
t86 = t106 * pkin(5) - t141 * t118 + t134;
t85 = t127 * t86 + t130 * t88;
t84 = -t127 * t88 + t130 * t86;
t1 = m(6) * (t89 ^ 2 + t90 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t87 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t124 ^ 2 * t132 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + t132 * Ifges(2,3) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + (-t87 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t116 * mrSges(3,1) - t117 * mrSges(3,2) + Ifges(3,3) * t122 / 0.2e1) * t122 + (t108 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,1) * t115 / 0.2e1) * t115 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t103 / 0.2e1) * t103 + ((-t116 * mrSges(3,3) + Ifges(3,5) * t122 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t129 / 0.2e1) * t138) * t129 + (-t101 * mrSges(4,1) + t102 * mrSges(4,2) + t117 * mrSges(3,3) - Ifges(4,5) * t115 + Ifges(3,6) * t122 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t129 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t131) * t138) * t131) * t138 + (-Ifges(4,6) * t135 - t108 * mrSges(4,1) + t102 * mrSges(4,3) + Ifges(4,4) * t115 + Ifges(4,2) * t114 / 0.2e1) * t114 + (t87 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t103 + Ifges(7,1) * t100 / 0.2e1) * t100 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + t89 * mrSges(6,2) - t90 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t118) * t118 + (t89 * mrSges(6,1) + t104 * mrSges(5,2) - t91 * mrSges(5,3) - t93 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t106 + (-Ifges(6,4) + Ifges(5,5)) * t118) * t106 + (t104 * mrSges(5,1) + t90 * mrSges(6,1) - t93 * mrSges(6,2) - t92 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t105 + (Ifges(6,5) - Ifges(5,6)) * t118 + (-Ifges(5,4) - Ifges(6,6)) * t106) * t105;
T  = t1;
