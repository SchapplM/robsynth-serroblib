% Calculate kinetic energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:03
% EndTime: 2019-03-09 11:24:04
% DurationCPUTime: 0.51s
% Computational Cost: add. (802->113), mult. (1788->169), div. (0->0), fcn. (1320->10), ass. (0->49)
t138 = -pkin(2) - pkin(9);
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t127 = sin(qJ(2));
t122 = sin(pkin(6));
t137 = qJD(1) * t122;
t133 = t127 * t137;
t114 = qJD(4) + t133;
t130 = cos(qJ(2));
t132 = -qJ(3) * t127 - pkin(1);
t105 = (t130 * t138 + t132) * t137;
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t116 = pkin(8) * t133;
t124 = cos(pkin(6));
t136 = qJD(1) * t124;
t120 = qJD(2) + t136;
t99 = qJD(3) + t116 + t138 * t120 + (-pkin(1) * t124 * t130 + pkin(3) * t122 * t127) * qJD(1);
t96 = t129 * t105 + t126 * t99;
t90 = qJ(5) * t114 + t96;
t134 = t130 * t137;
t135 = pkin(1) * t136;
t113 = pkin(8) * t134 + t127 * t135;
t107 = -t120 * qJ(3) - t113;
t104 = pkin(3) * t134 - t107;
t110 = t120 * t126 + t129 * t134;
t111 = t120 * t129 - t126 * t134;
t97 = pkin(4) * t110 - qJ(5) * t111 + t104;
t86 = t121 * t97 + t123 * t90;
t85 = -t121 * t90 + t123 * t97;
t95 = -t126 * t105 + t129 * t99;
t112 = t130 * t135 - t116;
t89 = -pkin(4) * t114 + qJD(5) - t95;
t131 = qJD(1) ^ 2;
t128 = cos(qJ(6));
t125 = sin(qJ(6));
t109 = qJD(6) + t110;
t108 = (-pkin(2) * t130 + t132) * t137;
t106 = -pkin(2) * t120 + qJD(3) - t112;
t101 = t111 * t123 + t114 * t121;
t100 = -t111 * t121 + t114 * t123;
t92 = t100 * t125 + t101 * t128;
t91 = t100 * t128 - t101 * t125;
t87 = -pkin(5) * t100 + t89;
t84 = pkin(10) * t100 + t86;
t83 = pkin(5) * t110 - pkin(10) * t101 + t85;
t82 = t125 * t83 + t128 * t84;
t81 = -t125 * t84 + t128 * t83;
t1 = t131 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t122 ^ 2 * t131 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + (t87 * mrSges(7,2) - t81 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (t95 * mrSges(5,1) - t96 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (t89 * mrSges(6,2) - t85 * mrSges(6,3) + Ifges(6,1) * t101 / 0.2e1) * t101 + (-t87 * mrSges(7,1) + t82 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (t104 * mrSges(5,2) - t95 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t111 / 0.2e1) * t111 + (-t89 * mrSges(6,1) + t86 * mrSges(6,3) + Ifges(6,4) * t101 + Ifges(6,2) * t100 / 0.2e1) * t100 + (t81 * mrSges(7,1) - t82 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t109 / 0.2e1) * t109 + ((-t107 * mrSges(4,1) + t108 * mrSges(4,2) + t113 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t130) * t137) * t130 + (t106 * mrSges(4,1) - t112 * mrSges(3,3) - t108 * mrSges(4,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t127 + (Ifges(3,4) + Ifges(4,6)) * t130) * t137) * t127) * t137 + (t112 * mrSges(3,1) - t113 * mrSges(3,2) + t106 * mrSges(4,2) - t107 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t120 + ((-Ifges(4,5) + Ifges(3,6)) * t130 + (-Ifges(4,4) + Ifges(3,5)) * t127) * t137) * t120 + (t104 * mrSges(5,1) + t85 * mrSges(6,1) - t86 * mrSges(6,2) - t96 * mrSges(5,3) - Ifges(5,4) * t111 + Ifges(6,5) * t101 - Ifges(5,6) * t114 + Ifges(6,6) * t100 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t110) * t110;
T  = t1;
