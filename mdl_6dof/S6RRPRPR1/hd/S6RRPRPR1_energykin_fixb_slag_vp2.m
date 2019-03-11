% Calculate kinetic energy for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:49
% EndTime: 2019-03-09 10:07:50
% DurationCPUTime: 0.65s
% Computational Cost: add. (958->110), mult. (2267->169), div. (0->0), fcn. (1728->10), ass. (0->43)
t131 = qJD(1) * (pkin(7) + qJ(3));
t129 = pkin(7) * mrSges(3,3);
t128 = cos(qJ(4));
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t115 = qJD(2) + qJD(4);
t121 = sin(qJ(4));
t122 = sin(qJ(2));
t112 = qJD(2) * pkin(2) - t122 * t131;
t124 = cos(qJ(2));
t113 = t124 * t131;
t117 = sin(pkin(10));
t119 = cos(pkin(10));
t103 = t119 * t112 - t113 * t117;
t110 = (t117 * t124 + t119 * t122) * qJD(1);
t96 = qJD(2) * pkin(3) - pkin(8) * t110 + t103;
t104 = t117 * t112 + t119 * t113;
t109 = (-t117 * t122 + t119 * t124) * qJD(1);
t97 = pkin(8) * t109 + t104;
t87 = t121 * t96 + t128 * t97;
t85 = qJ(5) * t115 + t87;
t101 = -t128 * t109 + t110 * t121;
t102 = t121 * t109 + t128 * t110;
t114 = qJD(3) + (-pkin(2) * t124 - pkin(1)) * qJD(1);
t105 = -pkin(3) * t109 + t114;
t90 = pkin(4) * t101 - qJ(5) * t102 + t105;
t81 = t116 * t90 + t118 * t85;
t80 = -t116 * t85 + t118 * t90;
t86 = -t121 * t97 + t128 * t96;
t84 = -t115 * pkin(4) + qJD(5) - t86;
t123 = cos(qJ(6));
t120 = sin(qJ(6));
t100 = qJD(6) + t101;
t99 = t102 * t118 + t115 * t116;
t98 = -t102 * t116 + t115 * t118;
t92 = t120 * t98 + t123 * t99;
t91 = -t120 * t99 + t123 * t98;
t82 = -t98 * pkin(5) + t84;
t79 = pkin(9) * t98 + t81;
t78 = pkin(5) * t101 - pkin(9) * t99 + t80;
t77 = t120 * t78 + t123 * t79;
t76 = -t120 * t79 + t123 * t78;
t1 = m(5) * (t105 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + (t84 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t99 / 0.2e1) * t99 + (t82 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,3) * t115 / 0.2e1) * t115 + (t114 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,1) * t110 / 0.2e1) * t110 + (-t84 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,2) * t98 / 0.2e1) * t98 + (-t82 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (-t114 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t110 + Ifges(4,2) * t109 / 0.2e1) * t109 + (t105 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,5) * t115 + Ifges(5,1) * t102 / 0.2e1) * t102 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,5) * t110 + Ifges(4,6) * t109 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t122 + Ifges(3,6) * t124 + (-mrSges(3,1) * t122 - mrSges(3,2) * t124) * pkin(7)) * qJD(1)) * qJD(2) + (t105 * mrSges(5,1) + t80 * mrSges(6,1) - t81 * mrSges(6,2) - t87 * mrSges(5,3) - Ifges(5,4) * t102 + Ifges(6,5) * t99 - Ifges(5,6) * t115 + Ifges(6,6) * t98 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t101) * t101 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t122 ^ 2 + t124 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t129) * t124) * t124 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t124 + (Ifges(3,1) / 0.2e1 + t129) * t122) * t122) * qJD(1) ^ 2;
T  = t1;
