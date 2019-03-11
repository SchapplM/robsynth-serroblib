% Calculate kinetic energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:05
% EndTime: 2019-03-09 09:01:06
% DurationCPUTime: 0.47s
% Computational Cost: add. (712->112), mult. (1924->165), div. (0->0), fcn. (1474->10), ass. (0->49)
t137 = pkin(3) + pkin(9);
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t119 = sin(pkin(11));
t124 = sin(qJ(2));
t127 = cos(qJ(2));
t120 = sin(pkin(6));
t135 = qJD(1) * t120;
t136 = cos(pkin(11));
t109 = (t119 * t127 + t124 * t136) * t135;
t134 = cos(pkin(6)) * qJD(1);
t118 = qJD(2) + t134;
t133 = pkin(1) * t134;
t117 = t127 * t133;
t132 = t124 * t135;
t103 = pkin(2) * t118 + t117 + (-pkin(8) - qJ(3)) * t132;
t131 = t127 * t135;
t111 = pkin(8) * t131 + t124 * t133;
t106 = qJ(3) * t131 + t111;
t94 = t103 * t136 - t119 * t106;
t130 = qJD(4) - t94;
t88 = t109 * pkin(4) - t118 * t137 + t130;
t108 = t119 * t132 - t131 * t136;
t112 = qJD(3) + (-pkin(2) * t127 - pkin(1)) * t135;
t129 = -qJ(4) * t109 + t112;
t91 = t108 * t137 + t129;
t85 = t123 * t88 + t126 * t91;
t95 = t119 * t103 + t136 * t106;
t93 = -t118 * qJ(4) - t95;
t84 = -t123 * t91 + t126 * t88;
t101 = t108 * t126 - t118 * t123;
t89 = -pkin(4) * t108 - t93;
t128 = qJD(1) ^ 2;
t125 = cos(qJ(6));
t122 = sin(qJ(6));
t110 = -pkin(8) * t132 + t117;
t107 = qJD(5) + t109;
t102 = t108 * t123 + t118 * t126;
t100 = qJD(6) - t101;
t98 = pkin(3) * t108 + t129;
t97 = t102 * t125 + t107 * t122;
t96 = -t102 * t122 + t107 * t125;
t92 = -t118 * pkin(3) + t130;
t86 = -pkin(5) * t101 - pkin(10) * t102 + t89;
t83 = pkin(10) * t107 + t85;
t82 = -pkin(5) * t107 - t84;
t81 = t122 * t86 + t125 * t83;
t80 = -t122 * t83 + t125 * t86;
t1 = m(3) * (pkin(1) ^ 2 * t120 ^ 2 * t128 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t128 * Ifges(2,3) / 0.2e1 + m(4) * (t112 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + (t82 * mrSges(7,2) - t80 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t84 * mrSges(6,1) - t85 * mrSges(6,2) + Ifges(6,3) * t107 / 0.2e1) * t107 + (-t82 * mrSges(7,1) + t81 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (t89 * mrSges(6,2) - t84 * mrSges(6,3) + Ifges(6,5) * t107 + Ifges(6,1) * t102 / 0.2e1) * t102 + (-t89 * mrSges(6,1) + t85 * mrSges(6,3) + Ifges(6,4) * t102 + Ifges(6,6) * t107 + Ifges(6,2) * t101 / 0.2e1) * t101 + (t80 * mrSges(7,1) - t81 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t92 * mrSges(5,1) + t112 * mrSges(4,2) - t94 * mrSges(4,3) - t98 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t109) * t109 + (t112 * mrSges(4,1) + t93 * mrSges(5,1) - t98 * mrSges(5,2) - t95 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t108 + (-Ifges(4,4) - Ifges(5,6)) * t109) * t108 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t127 / 0.2e1) * t127 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t127 + Ifges(3,1) * t124 / 0.2e1) * t124) * t135 + (-t110 * t124 + t111 * t127) * mrSges(3,3)) * t135 + (t110 * mrSges(3,1) + t94 * mrSges(4,1) - t111 * mrSges(3,2) - t95 * mrSges(4,2) + t92 * mrSges(5,2) - t93 * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t118 + (-Ifges(5,4) + Ifges(4,5)) * t109 + (Ifges(5,5) - Ifges(4,6)) * t108 + (Ifges(3,5) * t124 + Ifges(3,6) * t127) * t135) * t118;
T  = t1;
