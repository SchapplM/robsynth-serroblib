% Calculate kinetic energy for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:19
% EndTime: 2019-03-09 18:15:19
% DurationCPUTime: 0.63s
% Computational Cost: add. (1084->110), mult. (2231->170), div. (0->0), fcn. (1678->10), ass. (0->46)
t133 = -pkin(8) - pkin(7);
t132 = pkin(7) * mrSges(3,3);
t121 = sin(qJ(5));
t125 = cos(qJ(5));
t122 = sin(qJ(3));
t123 = sin(qJ(2));
t126 = cos(qJ(3));
t127 = cos(qJ(2));
t110 = (t122 * t127 + t123 * t126) * qJD(1);
t117 = qJD(2) + qJD(3);
t118 = sin(pkin(11));
t119 = cos(pkin(11));
t105 = t110 * t119 + t117 * t118;
t130 = qJD(1) * t127;
t131 = qJD(1) * t123;
t109 = t122 * t131 - t126 * t130;
t113 = qJD(2) * pkin(2) + t133 * t131;
t114 = t133 * t130;
t103 = t122 * t113 - t126 * t114;
t101 = qJ(4) * t117 + t103;
t115 = (-pkin(2) * t127 - pkin(1)) * qJD(1);
t98 = pkin(3) * t109 - qJ(4) * t110 + t115;
t91 = -t101 * t118 + t119 * t98;
t85 = pkin(4) * t109 - pkin(9) * t105 + t91;
t104 = -t110 * t118 + t117 * t119;
t92 = t119 * t101 + t118 * t98;
t90 = pkin(9) * t104 + t92;
t82 = t121 * t85 + t125 * t90;
t81 = -t121 * t90 + t125 * t85;
t102 = t113 * t126 + t122 * t114;
t108 = qJD(5) + t109;
t100 = -pkin(3) * t117 + qJD(4) - t102;
t93 = -pkin(4) * t104 + t100;
t124 = cos(qJ(6));
t120 = sin(qJ(6));
t106 = qJD(6) + t108;
t95 = t104 * t121 + t105 * t125;
t94 = t104 * t125 - t105 * t121;
t88 = t120 * t94 + t124 * t95;
t87 = -t120 * t95 + t124 * t94;
t86 = -pkin(5) * t94 + t93;
t80 = pkin(10) * t94 + t82;
t79 = pkin(5) * t108 - pkin(10) * t95 + t81;
t78 = t120 * t79 + t124 * t80;
t77 = -t120 * t80 + t124 * t79;
t1 = m(5) * (t100 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t93 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t115 ^ 2) / 0.2e1 + (t93 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,1) * t95 / 0.2e1) * t95 + (t86 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t88 / 0.2e1) * t88 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,3) * t117 / 0.2e1) * t117 + (t100 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (-t93 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t95 + Ifges(6,2) * t94 / 0.2e1) * t94 + (-t86 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t115 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,5) * t117 + Ifges(4,1) * t110 / 0.2e1) * t110 + (-t100 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,5) * t95 + Ifges(6,6) * t94 + Ifges(6,3) * t108 / 0.2e1) * t108 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t106 / 0.2e1) * t106 + (t115 * mrSges(4,1) + t91 * mrSges(5,1) - t92 * mrSges(5,2) - t103 * mrSges(4,3) - Ifges(4,4) * t110 + Ifges(5,5) * t105 - Ifges(4,6) * t117 + Ifges(5,6) * t104 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t109) * t109 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t123 ^ 2 + t127 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t132) * t127) * t127 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t127 + (Ifges(3,1) / 0.2e1 + t132) * t123) * t123) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t127 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t123) * qJD(2)) * qJD(1);
T  = t1;
