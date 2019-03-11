% Calculate kinetic energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:40
% EndTime: 2019-03-09 15:47:41
% DurationCPUTime: 0.49s
% Computational Cost: add. (932->111), mult. (2134->164), div. (0->0), fcn. (1666->10), ass. (0->48)
t138 = pkin(4) + pkin(10);
t120 = sin(pkin(11));
t137 = cos(pkin(11));
t135 = cos(pkin(6)) * qJD(1);
t119 = qJD(2) + t135;
t124 = sin(qJ(3));
t127 = cos(qJ(3));
t125 = sin(qJ(2));
t121 = sin(pkin(6));
t136 = qJD(1) * t121;
t133 = t125 * t136;
t112 = t119 * t124 + t127 * t133;
t128 = cos(qJ(2));
t132 = t128 * t136;
t115 = qJD(3) - t132;
t134 = pkin(1) * t135;
t114 = pkin(8) * t132 + t125 * t134;
t108 = pkin(9) * t119 + t114;
t110 = (-pkin(2) * t128 - pkin(9) * t125 - pkin(1)) * t136;
t98 = -t108 * t124 + t127 * t110;
t92 = pkin(3) * t115 - qJ(4) * t112 + t98;
t111 = t119 * t127 - t124 * t133;
t99 = t127 * t108 + t124 * t110;
t95 = qJ(4) * t111 + t99;
t89 = t120 * t92 + t137 * t95;
t87 = -qJ(5) * t115 - t89;
t88 = -t120 * t95 + t137 * t92;
t113 = -pkin(8) * t133 + t128 * t134;
t131 = qJD(5) - t88;
t103 = t120 * t111 + t137 * t112;
t107 = -pkin(2) * t119 - t113;
t101 = -pkin(3) * t111 + qJD(4) + t107;
t130 = -qJ(5) * t103 + t101;
t129 = qJD(1) ^ 2;
t126 = cos(qJ(6));
t123 = sin(qJ(6));
t102 = -t137 * t111 + t112 * t120;
t100 = qJD(6) + t103;
t97 = t102 * t123 + t115 * t126;
t96 = t102 * t126 - t115 * t123;
t90 = pkin(4) * t102 + t130;
t86 = -t115 * pkin(4) + t131;
t85 = t138 * t102 + t130;
t84 = -pkin(5) * t102 - t87;
t83 = t103 * pkin(5) - t138 * t115 + t131;
t82 = t123 * t83 + t126 * t85;
t81 = -t123 * t85 + t126 * t83;
t1 = m(4) * (t107 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t121 ^ 2 * t129 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t129 * Ifges(2,3) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t84 ^ 2) / 0.2e1 + (t84 * mrSges(7,2) - t81 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t107 * mrSges(4,2) - t98 * mrSges(4,3) + Ifges(4,1) * t112 / 0.2e1) * t112 + (-t84 * mrSges(7,1) + t82 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (-t107 * mrSges(4,1) + t99 * mrSges(4,3) + Ifges(4,4) * t112 + Ifges(4,2) * t111 / 0.2e1) * t111 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t128 / 0.2e1) * t128 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t128 + Ifges(3,1) * t125 / 0.2e1) * t125) * t136 + (-t113 * t125 + t114 * t128) * mrSges(3,3)) * t136 + (t113 * mrSges(3,1) - t114 * mrSges(3,2) + Ifges(3,3) * t119 / 0.2e1 + (Ifges(3,5) * t125 + Ifges(3,6) * t128) * t136) * t119 + (t81 * mrSges(7,1) - t82 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t86 * mrSges(6,1) + t101 * mrSges(5,2) - t88 * mrSges(5,3) - t90 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t103) * t103 + (t101 * mrSges(5,1) + t87 * mrSges(6,1) - t90 * mrSges(6,2) - t89 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102 + (-Ifges(5,4) - Ifges(6,6)) * t103) * t102 + (t98 * mrSges(4,1) + t88 * mrSges(5,1) - t99 * mrSges(4,2) - t89 * mrSges(5,2) + t86 * mrSges(6,2) - t87 * mrSges(6,3) + Ifges(4,5) * t112 + Ifges(4,6) * t111 + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t115 + (-Ifges(6,4) + Ifges(5,5)) * t103 + (Ifges(6,5) - Ifges(5,6)) * t102) * t115;
T  = t1;
