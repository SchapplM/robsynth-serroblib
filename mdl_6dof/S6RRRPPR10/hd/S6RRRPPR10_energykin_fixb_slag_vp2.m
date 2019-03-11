% Calculate kinetic energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:23
% EndTime: 2019-03-09 16:21:24
% DurationCPUTime: 0.49s
% Computational Cost: add. (866->111), mult. (1908->164), div. (0->0), fcn. (1454->10), ass. (0->48)
t140 = cos(qJ(3));
t139 = pkin(3) + qJ(5);
t122 = sin(pkin(11));
t124 = cos(pkin(11));
t138 = cos(pkin(6)) * qJD(1);
t121 = qJD(2) + t138;
t127 = sin(qJ(3));
t128 = sin(qJ(2));
t123 = sin(pkin(6));
t137 = t123 * qJD(1);
t135 = t128 * t137;
t112 = t127 * t121 + t140 * t135;
t130 = cos(qJ(2));
t134 = t130 * t137;
t117 = qJD(3) - t134;
t136 = pkin(1) * t138;
t114 = pkin(8) * t134 + t128 * t136;
t108 = pkin(9) * t121 + t114;
t109 = (-pkin(2) * t130 - pkin(9) * t128 - pkin(1)) * t137;
t100 = -t127 * t108 + t140 * t109;
t133 = qJD(4) - t100;
t91 = t112 * pkin(4) - t139 * t117 + t133;
t111 = -t140 * t121 + t127 * t135;
t113 = -pkin(8) * t135 + t130 * t136;
t107 = -pkin(2) * t121 - t113;
t132 = -qJ(4) * t112 + t107;
t93 = t139 * t111 + t132;
t87 = t122 * t91 + t124 * t93;
t101 = t140 * t108 + t127 * t109;
t99 = -t117 * qJ(4) - t101;
t86 = -t122 * t93 + t124 * t91;
t94 = -pkin(4) * t111 + qJD(5) - t99;
t131 = qJD(1) ^ 2;
t129 = cos(qJ(6));
t126 = sin(qJ(6));
t110 = qJD(6) + t112;
t103 = t111 * t122 + t117 * t124;
t102 = t111 * t124 - t117 * t122;
t98 = -t117 * pkin(3) + t133;
t97 = pkin(3) * t111 + t132;
t96 = t102 * t126 + t103 * t129;
t95 = t102 * t129 - t103 * t126;
t88 = -pkin(5) * t102 + t94;
t85 = pkin(10) * t102 + t87;
t84 = pkin(5) * t112 - pkin(10) * t103 + t86;
t83 = t126 * t84 + t129 * t85;
t82 = -t126 * t85 + t129 * t84;
t1 = m(3) * (pkin(1) ^ 2 * t123 ^ 2 * t131 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t131 * Ifges(2,3) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) / 0.2e1 + (t88 * mrSges(7,2) - t82 * mrSges(7,3) + Ifges(7,1) * t96 / 0.2e1) * t96 + (t94 * mrSges(6,2) - t86 * mrSges(6,3) + Ifges(6,1) * t103 / 0.2e1) * t103 + (-t88 * mrSges(7,1) + t83 * mrSges(7,3) + Ifges(7,4) * t96 + Ifges(7,2) * t95 / 0.2e1) * t95 + (-t94 * mrSges(6,1) + t87 * mrSges(6,3) + Ifges(6,4) * t103 + Ifges(6,2) * t102 / 0.2e1) * t102 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t130 / 0.2e1) * t130 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t130 + Ifges(3,1) * t128 / 0.2e1) * t128) * t137 + (-t113 * t128 + t114 * t130) * mrSges(3,3)) * t137 + (t113 * mrSges(3,1) - t114 * mrSges(3,2) + Ifges(3,3) * t121 / 0.2e1 + (Ifges(3,5) * t128 + Ifges(3,6) * t130) * t137) * t121 + (t82 * mrSges(7,1) - t83 * mrSges(7,2) + Ifges(7,5) * t96 + Ifges(7,6) * t95 + Ifges(7,3) * t110 / 0.2e1) * t110 + (t100 * mrSges(4,1) - t101 * mrSges(4,2) + t98 * mrSges(5,2) - t99 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t117) * t117 + (t107 * mrSges(4,1) + t99 * mrSges(5,1) - t97 * mrSges(5,2) - t101 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t111 + (Ifges(5,5) - Ifges(4,6)) * t117) * t111 + (t98 * mrSges(5,1) + t86 * mrSges(6,1) + t107 * mrSges(4,2) - t87 * mrSges(6,2) - t100 * mrSges(4,3) - t97 * mrSges(5,3) + Ifges(6,5) * t103 + Ifges(6,6) * t102 + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t112 + (-Ifges(5,4) + Ifges(4,5)) * t117 + (-Ifges(4,4) - Ifges(5,6)) * t111) * t112;
T  = t1;
