% Calculate kinetic energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S7RRRRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp2: Ifges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:13
% EndTime: 2019-03-10 06:29:14
% DurationCPUTime: 0.86s
% Computational Cost: add. (1046->103), mult. (2018->170), div. (0->0), fcn. (1684->12), ass. (0->47)
t118 = sin(qJ(3));
t124 = cos(qJ(3));
t119 = sin(qJ(2));
t128 = qJD(1) * t119;
t111 = -qJD(2) * t118 + t124 * t128;
t116 = sin(qJ(5));
t122 = cos(qJ(5));
t125 = cos(qJ(2));
t112 = qJD(1) * t125 + qJD(3);
t117 = sin(qJ(4));
t123 = cos(qJ(4));
t102 = t111 * t123 - t112 * t117;
t107 = t111 * pkin(2);
t97 = -pkin(3) * t102 - t107;
t110 = -qJD(2) * t124 - t118 * t128;
t108 = t110 * pkin(2);
t109 = qJD(4) + t110;
t98 = pkin(3) * t109 + t108 * t123;
t88 = t116 * t98 - t122 * t97;
t130 = t88 ^ 2;
t115 = sin(qJ(6));
t121 = cos(qJ(6));
t129 = t108 * t117;
t90 = t116 * t97 + t122 * t98;
t85 = t115 * t129 + t121 * t90;
t101 = -t111 * t117 - t112 * t123;
t100 = qJD(5) - t101;
t95 = t102 * t122 + t109 * t116;
t91 = t121 * t100 - t115 * t95;
t94 = -t102 * t116 + t109 * t122;
t120 = cos(qJ(7));
t114 = sin(qJ(7));
t106 = t108 ^ 2;
t105 = t107 ^ 2;
t104 = t117 ^ 2 * t106;
t93 = qJD(6) - t94;
t92 = t100 * t115 + t121 * t95;
t87 = qJD(7) + t91;
t84 = -t115 * t90 + t121 * t129;
t83 = t84 ^ 2;
t82 = -t114 * t93 + t120 * t92;
t81 = -t114 * t92 - t120 * t93;
t80 = -pkin(4) * t93 + t85;
t79 = pkin(4) * t92 + t88;
t78 = -t114 * t79 + t120 * t80;
t77 = -t114 * t80 - t120 * t79;
t1 = m(6) * (t90 ^ 2 + t104 + t130) / 0.2e1 + m(4) * (t106 + t105) / 0.2e1 + m(8) * (t77 ^ 2 + t78 ^ 2 + t83) / 0.2e1 + m(7) * (t85 ^ 2 + t130 + t83) / 0.2e1 + m(5) * (t106 * t123 ^ 2 + t104 + t105) / 0.2e1 + Ifges(5,3) * t109 ^ 2 / 0.2e1 + (t88 * mrSges(6,3) + Ifges(6,1) * t95 / 0.2e1) * t95 + (-t107 * mrSges(4,1) + Ifges(4,3) * t112 / 0.2e1) * t112 + (t107 * mrSges(4,3) + Ifges(4,5) * t112 + Ifges(4,1) * t111 / 0.2e1) * t111 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,3) * t93 / 0.2e1) * t93 + (t77 * mrSges(8,1) - t78 * mrSges(8,2) + Ifges(8,3) * t87 / 0.2e1) * t87 + (t90 * mrSges(6,3) + Ifges(6,4) * t95 + Ifges(6,2) * t94 / 0.2e1) * t94 + (t88 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,5) * t93 + Ifges(7,1) * t92 / 0.2e1) * t92 + (t84 * mrSges(8,2) - t77 * mrSges(8,3) + Ifges(8,5) * t87 + Ifges(8,1) * t82 / 0.2e1) * t82 + (-t107 * mrSges(5,2) + Ifges(5,5) * t109 + Ifges(5,1) * t102 / 0.2e1) * t102 + (-t88 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,6) * t93 + Ifges(7,2) * t91 / 0.2e1) * t91 + (-t84 * mrSges(8,1) + t78 * mrSges(8,3) + Ifges(8,4) * t82 + Ifges(8,6) * t87 + Ifges(8,2) * t81 / 0.2e1) * t81 + (Ifges(4,4) * t111 + Ifges(4,6) * t112 + Ifges(4,2) * t110 / 0.2e1) * t110 + (t107 * mrSges(5,1) + Ifges(5,4) * t102 + Ifges(5,6) * t109 + Ifges(5,2) * t101 / 0.2e1) * t101 + (-t112 * mrSges(4,2) + t110 * mrSges(4,3) + t123 * (-mrSges(5,2) * t109 + mrSges(5,3) * t101) + (-mrSges(5,1) * t109 - mrSges(6,1) * t94 + mrSges(6,2) * t95 + mrSges(5,3) * t102) * t117) * t108 + (-t88 * mrSges(6,1) - t90 * mrSges(6,2) + Ifges(6,5) * t95 + Ifges(6,6) * t94 + Ifges(6,3) * t100 / 0.2e1) * t100 + (Ifges(2,3) / 0.2e1 + Ifges(3,2) * t125 ^ 2 / 0.2e1 + (Ifges(3,4) * t125 + Ifges(3,1) * t119 / 0.2e1) * t119) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (Ifges(3,5) * t119 + Ifges(3,6) * t125) * qJD(1)) * qJD(2);
T  = t1;
