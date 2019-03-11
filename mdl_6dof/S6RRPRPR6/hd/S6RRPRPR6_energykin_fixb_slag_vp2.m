% Calculate kinetic energy for
% S6RRPRPR6
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:02
% EndTime: 2019-03-09 10:38:03
% DurationCPUTime: 0.59s
% Computational Cost: add. (778->112), mult. (2078->166), div. (0->0), fcn. (1618->10), ass. (0->49)
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t120 = sin(pkin(6));
t136 = qJD(1) * t120;
t109 = (t119 * t125 - t121 * t127) * t136;
t138 = pkin(4) + pkin(10);
t137 = cos(qJ(4));
t124 = sin(qJ(4));
t135 = cos(pkin(6)) * qJD(1);
t118 = qJD(2) + t135;
t134 = pkin(1) * t135;
t117 = t127 * t134;
t133 = t125 * t136;
t104 = pkin(2) * t118 + t117 + (-pkin(8) - qJ(3)) * t133;
t132 = t127 * t136;
t112 = pkin(8) * t132 + t125 * t134;
t107 = qJ(3) * t132 + t112;
t95 = t119 * t104 + t121 * t107;
t93 = pkin(9) * t118 + t95;
t110 = (t119 * t127 + t121 * t125) * t136;
t113 = qJD(3) + (-pkin(2) * t127 - pkin(1)) * t136;
t99 = pkin(3) * t109 - pkin(9) * t110 + t113;
t89 = t124 * t99 + t137 * t93;
t94 = t104 * t121 - t119 * t107;
t108 = qJD(4) + t109;
t86 = -qJ(5) * t108 - t89;
t88 = -t124 * t93 + t137 * t99;
t130 = qJD(5) - t88;
t92 = -pkin(3) * t118 - t94;
t103 = t137 * t110 + t124 * t118;
t129 = -qJ(5) * t103 + t92;
t128 = qJD(1) ^ 2;
t126 = cos(qJ(6));
t123 = sin(qJ(6));
t111 = -pkin(8) * t133 + t117;
t102 = t110 * t124 - t137 * t118;
t101 = qJD(6) + t103;
t97 = t102 * t123 + t108 * t126;
t96 = t102 * t126 - t108 * t123;
t87 = pkin(4) * t102 + t129;
t85 = -t108 * pkin(4) + t130;
t84 = t138 * t102 + t129;
t83 = -pkin(5) * t102 - t86;
t82 = t103 * pkin(5) - t138 * t108 + t130;
t81 = t123 * t82 + t126 * t84;
t80 = -t123 * t84 + t126 * t82;
t1 = m(3) * (pkin(1) ^ 2 * t120 ^ 2 * t128 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t128 * Ifges(2,3) / 0.2e1 + m(4) * (t113 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + (t83 * mrSges(7,2) - t80 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t113 * mrSges(4,2) - t94 * mrSges(4,3) + Ifges(4,1) * t110 / 0.2e1) * t110 + (-t83 * mrSges(7,1) + t81 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 - (-t113 * mrSges(4,1) + t95 * mrSges(4,3) + Ifges(4,4) * t110 - Ifges(4,2) * t109 / 0.2e1) * t109 + (t80 * mrSges(7,1) - t81 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t101 / 0.2e1) * t101 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + t85 * mrSges(6,2) - t86 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t108) * t108 + (t85 * mrSges(6,1) + t92 * mrSges(5,2) - t88 * mrSges(5,3) - t87 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t103 + (-Ifges(6,4) + Ifges(5,5)) * t108) * t103 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t127 / 0.2e1) * t127 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t127 + Ifges(3,1) * t125 / 0.2e1) * t125) * t136 + (-t111 * t125 + t112 * t127) * mrSges(3,3)) * t136 + (t111 * mrSges(3,1) + t94 * mrSges(4,1) - t112 * mrSges(3,2) - t95 * mrSges(4,2) + Ifges(4,5) * t110 - Ifges(4,6) * t109 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t118 + (Ifges(3,5) * t125 + Ifges(3,6) * t127) * t136) * t118 + (t92 * mrSges(5,1) + t86 * mrSges(6,1) - t87 * mrSges(6,2) - t89 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102 + (Ifges(6,5) - Ifges(5,6)) * t108 + (-Ifges(5,4) - Ifges(6,6)) * t103) * t102;
T  = t1;
