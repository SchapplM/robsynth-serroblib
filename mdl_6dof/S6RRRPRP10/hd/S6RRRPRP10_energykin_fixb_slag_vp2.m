% Calculate kinetic energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:23
% EndTime: 2019-03-09 17:29:23
% DurationCPUTime: 0.53s
% Computational Cost: add. (1076->110), mult. (2404->164), div. (0->0), fcn. (1896->10), ass. (0->45)
t138 = cos(qJ(5));
t127 = sin(qJ(5));
t137 = cos(pkin(6)) * qJD(1);
t122 = qJD(2) + t137;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t129 = sin(qJ(2));
t124 = sin(pkin(6));
t136 = t124 * qJD(1);
t133 = t129 * t136;
t114 = t122 * t128 + t130 * t133;
t131 = cos(qJ(2));
t134 = t131 * t136;
t118 = qJD(3) - t134;
t123 = sin(pkin(11));
t125 = cos(pkin(11));
t105 = t114 * t125 + t118 * t123;
t113 = -t130 * t122 + t128 * t133;
t135 = pkin(1) * t137;
t116 = pkin(8) * t134 + t129 * t135;
t110 = pkin(9) * t122 + t116;
t111 = (-pkin(2) * t131 - pkin(9) * t129 - pkin(1)) * t136;
t102 = t130 * t110 + t128 * t111;
t100 = qJ(4) * t118 + t102;
t115 = -pkin(8) * t133 + t131 * t135;
t109 = -pkin(2) * t122 - t115;
t97 = pkin(3) * t113 - qJ(4) * t114 + t109;
t90 = -t100 * t123 + t125 * t97;
t87 = pkin(4) * t113 - pkin(10) * t105 + t90;
t104 = -t114 * t123 + t118 * t125;
t91 = t125 * t100 + t123 * t97;
t89 = pkin(10) * t104 + t91;
t84 = t127 * t87 + t138 * t89;
t101 = -t128 * t110 + t111 * t130;
t83 = -t127 * t89 + t138 * t87;
t99 = -pkin(3) * t118 + qJD(4) - t101;
t92 = -pkin(4) * t104 + t99;
t132 = qJD(1) ^ 2;
t112 = qJD(5) + t113;
t94 = t127 * t104 + t105 * t138;
t93 = -t104 * t138 + t105 * t127;
t85 = pkin(5) * t93 - qJ(6) * t94 + t92;
t82 = qJ(6) * t112 + t84;
t81 = -t112 * pkin(5) + qJD(6) - t83;
t1 = m(7) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t124 ^ 2 * t132 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + t132 * Ifges(2,3) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t92 ^ 2) / 0.2e1 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * t118 / 0.2e1) * t118 + (t99 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (t109 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,5) * t118 + Ifges(4,1) * t114 / 0.2e1) * t114 + (-t99 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t131 / 0.2e1) * t131 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t131 + Ifges(3,1) * t129 / 0.2e1) * t129) * t136 + (-t115 * t129 + t116 * t131) * mrSges(3,3)) * t136 + (t115 * mrSges(3,1) - t116 * mrSges(3,2) + Ifges(3,3) * t122 / 0.2e1 + (Ifges(3,5) * t129 + Ifges(3,6) * t131) * t136) * t122 + (t92 * mrSges(6,2) + t81 * mrSges(7,2) - t83 * mrSges(6,3) - t85 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t94) * t94 + (t92 * mrSges(6,1) + t85 * mrSges(7,1) - t82 * mrSges(7,2) - t84 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t93 + (-Ifges(6,4) + Ifges(7,5)) * t94) * t93 + (t109 * mrSges(4,1) + t90 * mrSges(5,1) - t91 * mrSges(5,2) - t102 * mrSges(4,3) - Ifges(4,4) * t114 + Ifges(5,5) * t105 - Ifges(4,6) * t118 + Ifges(5,6) * t104 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t113) * t113 + (t83 * mrSges(6,1) - t81 * mrSges(7,1) - t84 * mrSges(6,2) + t82 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t112 + (Ifges(7,4) + Ifges(6,5)) * t94 + (-Ifges(6,6) + Ifges(7,6)) * t93) * t112;
T  = t1;
