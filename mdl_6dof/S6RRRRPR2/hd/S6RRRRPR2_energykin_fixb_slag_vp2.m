% Calculate kinetic energy for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:19
% EndTime: 2019-03-09 21:57:20
% DurationCPUTime: 0.65s
% Computational Cost: add. (1036->110), mult. (2267->170), div. (0->0), fcn. (1728->10), ass. (0->44)
t132 = qJD(1) * (-pkin(8) - pkin(7));
t130 = pkin(7) * mrSges(3,3);
t129 = cos(qJ(4));
t118 = sin(pkin(11));
t119 = cos(pkin(11));
t117 = qJD(2) + qJD(3);
t116 = qJD(4) + t117;
t121 = sin(qJ(4));
t123 = sin(qJ(2));
t113 = qJD(2) * pkin(2) + t123 * t132;
t126 = cos(qJ(2));
t114 = t126 * t132;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t104 = t125 * t113 + t114 * t122;
t111 = (t122 * t126 + t123 * t125) * qJD(1);
t95 = pkin(3) * t117 - pkin(9) * t111 + t104;
t105 = t122 * t113 - t125 * t114;
t110 = (-t122 * t123 + t125 * t126) * qJD(1);
t98 = pkin(9) * t110 + t105;
t88 = t121 * t95 + t129 * t98;
t86 = qJ(5) * t116 + t88;
t102 = -t129 * t110 + t111 * t121;
t103 = t121 * t110 + t129 * t111;
t115 = (-pkin(2) * t126 - pkin(1)) * qJD(1);
t106 = -pkin(3) * t110 + t115;
t93 = pkin(4) * t102 - qJ(5) * t103 + t106;
t82 = t118 * t93 + t119 * t86;
t81 = -t118 * t86 + t119 * t93;
t87 = -t121 * t98 + t129 * t95;
t85 = -t116 * pkin(4) + qJD(5) - t87;
t124 = cos(qJ(6));
t120 = sin(qJ(6));
t101 = qJD(6) + t102;
t100 = t103 * t119 + t116 * t118;
t99 = -t103 * t118 + t116 * t119;
t92 = t100 * t124 + t120 * t99;
t91 = -t100 * t120 + t124 * t99;
t83 = -t99 * pkin(5) + t85;
t80 = pkin(10) * t99 + t82;
t79 = pkin(5) * t102 - pkin(10) * t100 + t81;
t78 = t120 * t79 + t124 * t80;
t77 = -t120 * t80 + t124 * t79;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t104 ^ 2 + t105 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t83 ^ 2) / 0.2e1 + (-t85 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,2) * t99 / 0.2e1) * t99 + (t83 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (t104 * mrSges(4,1) - t105 * mrSges(4,2) + Ifges(4,3) * t117 / 0.2e1) * t117 + (t87 * mrSges(5,1) - t88 * mrSges(5,2) + Ifges(5,3) * t116 / 0.2e1) * t116 + (-t83 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (t115 * mrSges(4,2) - t104 * mrSges(4,3) + Ifges(4,5) * t117 + Ifges(4,1) * t111 / 0.2e1) * t111 + (t106 * mrSges(5,2) - t87 * mrSges(5,3) + Ifges(5,5) * t116 + Ifges(5,1) * t103 / 0.2e1) * t103 + (t85 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,1) * t100 / 0.2e1) * t100 + (-t115 * mrSges(4,1) + t105 * mrSges(4,3) + Ifges(4,4) * t111 + Ifges(4,6) * t117 + Ifges(4,2) * t110 / 0.2e1) * t110 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t101 / 0.2e1) * t101 + (t106 * mrSges(5,1) + t81 * mrSges(6,1) - t82 * mrSges(6,2) - t88 * mrSges(5,3) - Ifges(5,4) * t103 + Ifges(6,5) * t100 - Ifges(5,6) * t116 + Ifges(6,6) * t99 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102) * t102 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t123 ^ 2 + t126 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t130) * t126) * t126 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t126 + (Ifges(3,1) / 0.2e1 + t130) * t123) * t123) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t126 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t123) * qJD(2)) * qJD(1);
T  = t1;
