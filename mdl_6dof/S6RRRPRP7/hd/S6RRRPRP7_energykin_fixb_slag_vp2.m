% Calculate kinetic energy for
% S6RRRPRP7
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:55
% EndTime: 2019-03-09 17:03:56
% DurationCPUTime: 0.50s
% Computational Cost: add. (1062->110), mult. (2420->164), div. (0->0), fcn. (1914->10), ass. (0->45)
t136 = cos(qJ(5));
t125 = sin(qJ(5));
t129 = cos(qJ(2));
t122 = sin(pkin(6));
t134 = t122 * qJD(1);
t131 = t129 * t134;
t116 = qJD(3) - t131;
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t135 = cos(pkin(6)) * qJD(1);
t120 = qJD(2) + t135;
t126 = sin(qJ(3));
t128 = cos(qJ(3));
t127 = sin(qJ(2));
t132 = t127 * t134;
t112 = t120 * t126 + t128 * t132;
t133 = pkin(1) * t135;
t114 = pkin(8) * t131 + t127 * t133;
t109 = pkin(9) * t120 + t114;
t110 = (-pkin(2) * t129 - pkin(9) * t127 - pkin(1)) * t134;
t99 = -t109 * t126 + t128 * t110;
t93 = pkin(3) * t116 - qJ(4) * t112 + t99;
t100 = t128 * t109 + t126 * t110;
t111 = t120 * t128 - t126 * t132;
t96 = qJ(4) * t111 + t100;
t89 = t121 * t93 + t123 * t96;
t87 = pkin(10) * t116 + t89;
t113 = -pkin(8) * t132 + t129 * t133;
t108 = -pkin(2) * t120 - t113;
t102 = -pkin(3) * t111 + qJD(4) + t108;
t103 = t111 * t123 - t112 * t121;
t104 = t111 * t121 + t112 * t123;
t91 = -pkin(4) * t103 - pkin(10) * t104 + t102;
t83 = t125 * t91 + t136 * t87;
t88 = -t121 * t96 + t123 * t93;
t86 = -pkin(4) * t116 - t88;
t82 = -t125 * t87 + t136 * t91;
t130 = qJD(1) ^ 2;
t101 = qJD(5) - t103;
t98 = t104 * t136 + t125 * t116;
t97 = t104 * t125 - t116 * t136;
t84 = pkin(5) * t97 - qJ(6) * t98 + t86;
t81 = qJ(6) * t101 + t83;
t80 = -t101 * pkin(5) + qJD(6) - t82;
t1 = m(7) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t122 ^ 2 * t130 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t130 * Ifges(2,3) / 0.2e1 + m(4) * (t100 ^ 2 + t108 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + (t108 * mrSges(4,2) - t99 * mrSges(4,3) + Ifges(4,1) * t112 / 0.2e1) * t112 + (t102 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,1) * t104 / 0.2e1) * t104 + (-t108 * mrSges(4,1) + t100 * mrSges(4,3) + Ifges(4,4) * t112 + Ifges(4,2) * t111 / 0.2e1) * t111 + (-t102 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,2) * t103 / 0.2e1) * t103 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t129 / 0.2e1) * t129 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t129 + Ifges(3,1) * t127 / 0.2e1) * t127) * t134 + (-t113 * t127 + t114 * t129) * mrSges(3,3)) * t134 + (t113 * mrSges(3,1) - t114 * mrSges(3,2) + Ifges(3,3) * t120 / 0.2e1 + (Ifges(3,5) * t127 + Ifges(3,6) * t129) * t134) * t120 + (t86 * mrSges(6,2) + t80 * mrSges(7,2) - t82 * mrSges(6,3) - t84 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t98) * t98 + (t86 * mrSges(6,1) + t84 * mrSges(7,1) - t81 * mrSges(7,2) - t83 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t97 + (-Ifges(6,4) + Ifges(7,5)) * t98) * t97 + (t99 * mrSges(4,1) + t88 * mrSges(5,1) - t100 * mrSges(4,2) - t89 * mrSges(5,2) + Ifges(4,5) * t112 + Ifges(5,5) * t104 + Ifges(4,6) * t111 + Ifges(5,6) * t103 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t116) * t116 + (t82 * mrSges(6,1) - t80 * mrSges(7,1) - t83 * mrSges(6,2) + t81 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t101 + (Ifges(7,4) + Ifges(6,5)) * t98 + (-Ifges(6,6) + Ifges(7,6)) * t97) * t101;
T  = t1;
