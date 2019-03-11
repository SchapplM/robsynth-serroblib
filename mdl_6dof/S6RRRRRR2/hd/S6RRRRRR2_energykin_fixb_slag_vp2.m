% Calculate kinetic energy for
% S6RRRRRR2
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:07
% EndTime: 2019-03-10 03:33:08
% DurationCPUTime: 0.71s
% Computational Cost: add. (1048->110), mult. (2267->172), div. (0->0), fcn. (1728->10), ass. (0->45)
t134 = qJD(1) * (-pkin(8) - pkin(7));
t132 = pkin(7) * mrSges(3,3);
t121 = sin(qJ(5));
t126 = cos(qJ(5));
t119 = qJD(2) + qJD(3);
t118 = qJD(4) + t119;
t122 = sin(qJ(4));
t127 = cos(qJ(4));
t124 = sin(qJ(2));
t115 = qJD(2) * pkin(2) + t124 * t134;
t129 = cos(qJ(2));
t116 = t129 * t134;
t123 = sin(qJ(3));
t128 = cos(qJ(3));
t106 = t128 * t115 + t116 * t123;
t113 = (t123 * t129 + t124 * t128) * qJD(1);
t96 = pkin(3) * t119 - pkin(9) * t113 + t106;
t107 = t123 * t115 - t128 * t116;
t112 = (-t123 * t124 + t128 * t129) * qJD(1);
t99 = pkin(9) * t112 + t107;
t89 = t122 * t96 + t127 * t99;
t87 = pkin(10) * t118 + t89;
t104 = t112 * t127 - t122 * t113;
t105 = t112 * t122 + t113 * t127;
t117 = (-pkin(2) * t129 - pkin(1)) * qJD(1);
t108 = -pkin(3) * t112 + t117;
t94 = -pkin(4) * t104 - pkin(10) * t105 + t108;
t83 = t121 * t94 + t126 * t87;
t82 = -t121 * t87 + t126 * t94;
t88 = -t122 * t99 + t127 * t96;
t103 = qJD(5) - t104;
t86 = -pkin(4) * t118 - t88;
t125 = cos(qJ(6));
t120 = sin(qJ(6));
t102 = qJD(6) + t103;
t101 = t105 * t126 + t118 * t121;
t100 = -t105 * t121 + t118 * t126;
t93 = t100 * t120 + t101 * t125;
t92 = t100 * t125 - t101 * t120;
t84 = -pkin(5) * t100 + t86;
t81 = pkin(11) * t100 + t83;
t80 = pkin(5) * t103 - pkin(11) * t101 + t82;
t79 = t120 * t80 + t125 * t81;
t78 = -t120 * t81 + t125 * t80;
t1 = m(4) * (t106 ^ 2 + t107 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t78 ^ 2 + t79 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t84 * mrSges(7,2) - t78 * mrSges(7,3) + Ifges(7,1) * t93 / 0.2e1) * t93 + (t106 * mrSges(4,1) - t107 * mrSges(4,2) + Ifges(4,3) * t119 / 0.2e1) * t119 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + Ifges(5,3) * t118 / 0.2e1) * t118 + (t82 * mrSges(6,1) - t83 * mrSges(6,2) + Ifges(6,3) * t103 / 0.2e1) * t103 + (-t84 * mrSges(7,1) + t79 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,2) * t92 / 0.2e1) * t92 + (t117 * mrSges(4,2) - t106 * mrSges(4,3) + Ifges(4,5) * t119 + Ifges(4,1) * t113 / 0.2e1) * t113 + (t108 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,5) * t118 + Ifges(5,1) * t105 / 0.2e1) * t105 + (t86 * mrSges(6,2) - t82 * mrSges(6,3) + Ifges(6,5) * t103 + Ifges(6,1) * t101 / 0.2e1) * t101 + (-t117 * mrSges(4,1) + t107 * mrSges(4,3) + Ifges(4,4) * t113 + Ifges(4,6) * t119 + Ifges(4,2) * t112 / 0.2e1) * t112 + (-t108 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,6) * t118 + Ifges(5,2) * t104 / 0.2e1) * t104 + (t78 * mrSges(7,1) - t79 * mrSges(7,2) + Ifges(7,5) * t93 + Ifges(7,6) * t92 + Ifges(7,3) * t102 / 0.2e1) * t102 + (-t86 * mrSges(6,1) + t83 * mrSges(6,3) + Ifges(6,4) * t101 + Ifges(6,6) * t103 + Ifges(6,2) * t100 / 0.2e1) * t100 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t124 ^ 2 + t129 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t132) * t129) * t129 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t129 + (Ifges(3,1) / 0.2e1 + t132) * t124) * t124) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t129 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t124) * qJD(2)) * qJD(1);
T  = t1;
