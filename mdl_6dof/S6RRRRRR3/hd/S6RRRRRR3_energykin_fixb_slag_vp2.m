% Calculate kinetic energy for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:07
% EndTime: 2019-03-10 03:39:08
% DurationCPUTime: 0.65s
% Computational Cost: add. (1112->110), mult. (2231->172), div. (0->0), fcn. (1678->10), ass. (0->47)
t135 = -pkin(8) - pkin(7);
t134 = pkin(7) * mrSges(3,3);
t121 = sin(qJ(5));
t126 = cos(qJ(5));
t123 = sin(qJ(3));
t124 = sin(qJ(2));
t128 = cos(qJ(3));
t129 = cos(qJ(2));
t112 = (t123 * t129 + t124 * t128) * qJD(1);
t119 = qJD(2) + qJD(3);
t122 = sin(qJ(4));
t127 = cos(qJ(4));
t106 = t112 * t127 + t119 * t122;
t132 = qJD(1) * t129;
t133 = qJD(1) * t124;
t111 = -t123 * t133 + t128 * t132;
t110 = qJD(4) - t111;
t115 = qJD(2) * pkin(2) + t135 * t133;
t116 = t135 * t132;
t104 = t123 * t115 - t128 * t116;
t102 = pkin(9) * t119 + t104;
t117 = (-pkin(2) * t129 - pkin(1)) * qJD(1);
t99 = -pkin(3) * t111 - pkin(9) * t112 + t117;
t92 = -t102 * t122 + t127 * t99;
t86 = pkin(4) * t110 - pkin(10) * t106 + t92;
t105 = -t112 * t122 + t119 * t127;
t93 = t127 * t102 + t122 * t99;
t91 = pkin(10) * t105 + t93;
t83 = t121 * t86 + t126 * t91;
t82 = -t121 * t91 + t126 * t86;
t103 = t115 * t128 + t123 * t116;
t101 = -pkin(3) * t119 - t103;
t108 = qJD(5) + t110;
t94 = -pkin(4) * t105 + t101;
t125 = cos(qJ(6));
t120 = sin(qJ(6));
t107 = qJD(6) + t108;
t96 = t105 * t121 + t106 * t126;
t95 = t105 * t126 - t106 * t121;
t89 = t120 * t95 + t125 * t96;
t88 = -t120 * t96 + t125 * t95;
t87 = -pkin(5) * t95 + t94;
t81 = pkin(11) * t95 + t83;
t80 = pkin(5) * t108 - pkin(11) * t96 + t82;
t79 = t120 * t80 + t125 * t81;
t78 = -t120 * t81 + t125 * t80;
t1 = m(4) * (t103 ^ 2 + t104 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t78 ^ 2 + t79 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t94 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(5) * (t101 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + (t94 * mrSges(6,2) - t82 * mrSges(6,3) + Ifges(6,1) * t96 / 0.2e1) * t96 + (t87 * mrSges(7,2) - t78 * mrSges(7,3) + Ifges(7,1) * t89 / 0.2e1) * t89 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,3) * t119 / 0.2e1) * t119 + (t92 * mrSges(5,1) - t93 * mrSges(5,2) + Ifges(5,3) * t110 / 0.2e1) * t110 + (-t94 * mrSges(6,1) + t83 * mrSges(6,3) + Ifges(6,4) * t96 + Ifges(6,2) * t95 / 0.2e1) * t95 + (-t87 * mrSges(7,1) + t79 * mrSges(7,3) + Ifges(7,4) * t89 + Ifges(7,2) * t88 / 0.2e1) * t88 + (t117 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,5) * t119 + Ifges(4,1) * t112 / 0.2e1) * t112 + (t101 * mrSges(5,2) - t92 * mrSges(5,3) + Ifges(5,5) * t110 + Ifges(5,1) * t106 / 0.2e1) * t106 + (-t117 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t112 + Ifges(4,6) * t119 + Ifges(4,2) * t111 / 0.2e1) * t111 + (t82 * mrSges(6,1) - t83 * mrSges(6,2) + Ifges(6,5) * t96 + Ifges(6,6) * t95 + Ifges(6,3) * t108 / 0.2e1) * t108 + (t78 * mrSges(7,1) - t79 * mrSges(7,2) + Ifges(7,5) * t89 + Ifges(7,6) * t88 + Ifges(7,3) * t107 / 0.2e1) * t107 + (-t101 * mrSges(5,1) + t93 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,6) * t110 + Ifges(5,2) * t105 / 0.2e1) * t105 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t124 ^ 2 + t129 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t134) * t129) * t129 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t129 + (Ifges(3,1) / 0.2e1 + t134) * t124) * t124) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t129 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t124) * qJD(2)) * qJD(1);
T  = t1;
