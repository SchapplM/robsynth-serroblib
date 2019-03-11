% Calculate kinetic energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:38:47
% EndTime: 2019-03-09 09:38:47
% DurationCPUTime: 0.50s
% Computational Cost: add. (800->113), mult. (1826->169), div. (0->0), fcn. (1350->10), ass. (0->49)
t135 = -pkin(2) - qJ(4);
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t121 = cos(pkin(6));
t133 = qJD(1) * t121;
t117 = qJD(2) + t133;
t118 = sin(pkin(11));
t120 = cos(pkin(11));
t127 = cos(qJ(2));
t119 = sin(pkin(6));
t134 = qJD(1) * t119;
t130 = t127 * t134;
t108 = t117 * t120 - t118 * t130;
t124 = sin(qJ(2));
t131 = t124 * t134;
t113 = pkin(8) * t131;
t100 = qJD(3) + t113 + t135 * t117 + (-pkin(1) * t121 * t127 + pkin(3) * t119 * t124) * qJD(1);
t129 = -qJ(3) * t124 - pkin(1);
t103 = (t127 * t135 + t129) * t134;
t90 = t120 * t100 - t103 * t118;
t87 = pkin(4) * t131 - pkin(9) * t108 + t90;
t107 = -t117 * t118 - t120 * t130;
t91 = t118 * t100 + t120 * t103;
t89 = pkin(9) * t107 + t91;
t84 = t123 * t87 + t126 * t89;
t132 = pkin(1) * t133;
t110 = pkin(8) * t130 + t124 * t132;
t105 = -t117 * qJ(3) - t110;
t83 = -t123 * t89 + t126 * t87;
t101 = pkin(3) * t130 + qJD(4) - t105;
t96 = t107 * t126 - t108 * t123;
t109 = t127 * t132 - t113;
t94 = -pkin(4) * t107 + t101;
t128 = qJD(1) ^ 2;
t125 = cos(qJ(6));
t122 = sin(qJ(6));
t111 = qJD(5) + t131;
t106 = (-pkin(2) * t127 + t129) * t134;
t104 = -pkin(2) * t117 + qJD(3) - t109;
t97 = t107 * t123 + t108 * t126;
t95 = qJD(6) - t96;
t93 = t111 * t122 + t125 * t97;
t92 = t111 * t125 - t122 * t97;
t85 = -pkin(5) * t96 - pkin(10) * t97 + t94;
t82 = pkin(10) * t111 + t84;
t81 = -pkin(5) * t111 - t83;
t80 = t122 * t85 + t125 * t82;
t79 = -t122 * t82 + t125 * t85;
t1 = m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t128 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t119 ^ 2 * t128 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + (t94 * mrSges(6,2) - t83 * mrSges(6,3) + Ifges(6,1) * t97 / 0.2e1) * t97 + (t79 * mrSges(7,1) - t80 * mrSges(7,2) + Ifges(7,3) * t95 / 0.2e1) * t95 + (t101 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,1) * t108 / 0.2e1) * t108 + (-t94 * mrSges(6,1) + t84 * mrSges(6,3) + Ifges(6,4) * t97 + Ifges(6,2) * t96 / 0.2e1) * t96 + (t81 * mrSges(7,2) - t79 * mrSges(7,3) + Ifges(7,5) * t95 + Ifges(7,1) * t93 / 0.2e1) * t93 + (-t101 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t108 + Ifges(5,2) * t107 / 0.2e1) * t107 + (-t81 * mrSges(7,1) + t80 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,6) * t95 + Ifges(7,2) * t92 / 0.2e1) * t92 + (t83 * mrSges(6,1) - t84 * mrSges(6,2) + Ifges(6,5) * t97 + Ifges(6,6) * t96 + Ifges(6,3) * t111 / 0.2e1) * t111 + ((-t105 * mrSges(4,1) + t106 * mrSges(4,2) + t110 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t127) * t134) * t127 + (t104 * mrSges(4,1) + t90 * mrSges(5,1) - t91 * mrSges(5,2) - t109 * mrSges(3,3) - t106 * mrSges(4,3) + Ifges(5,5) * t108 + Ifges(5,6) * t107 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t124 + (Ifges(3,4) + Ifges(4,6)) * t127) * t134) * t124) * t134 + (t109 * mrSges(3,1) - t110 * mrSges(3,2) + t104 * mrSges(4,2) - t105 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t117 + ((-Ifges(4,5) + Ifges(3,6)) * t127 + (-Ifges(4,4) + Ifges(3,5)) * t124) * t134) * t117;
T  = t1;
