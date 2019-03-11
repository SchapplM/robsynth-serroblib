% Calculate kinetic energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:36
% EndTime: 2019-03-09 13:55:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (598->108), mult. (1161->158), div. (0->0), fcn. (738->8), ass. (0->40)
t130 = pkin(7) * mrSges(3,3);
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t119 = sin(qJ(4));
t123 = cos(qJ(4));
t124 = cos(qJ(2));
t128 = qJD(1) * t124;
t120 = sin(qJ(2));
t129 = qJD(1) * t120;
t101 = -t119 * t129 - t123 * t128;
t102 = (-t119 * t124 + t120 * t123) * qJD(1);
t104 = -qJD(1) * pkin(1) - pkin(2) * t128 - qJ(3) * t129;
t95 = pkin(3) * t128 - t104;
t86 = -pkin(4) * t101 - pkin(9) * t102 + t95;
t114 = -qJD(2) + qJD(4);
t106 = pkin(7) * t128 + qJD(2) * qJ(3);
t103 = -pkin(8) * t128 + t106;
t127 = pkin(7) * t129 + qJD(3);
t96 = -pkin(8) * t129 + (-pkin(2) - pkin(3)) * qJD(2) + t127;
t91 = t123 * t103 + t119 * t96;
t89 = pkin(9) * t114 + t91;
t80 = t118 * t86 + t122 * t89;
t100 = qJD(5) - t101;
t79 = -t118 * t89 + t122 * t86;
t90 = -t119 * t103 + t123 * t96;
t88 = -pkin(4) * t114 - t90;
t121 = cos(qJ(6));
t117 = sin(qJ(6));
t105 = -qJD(2) * pkin(2) + t127;
t99 = qJD(6) + t100;
t93 = t102 * t122 + t114 * t118;
t92 = -t102 * t118 + t114 * t122;
t83 = t117 * t92 + t121 * t93;
t82 = -t117 * t93 + t121 * t92;
t81 = -pkin(5) * t92 + t88;
t78 = pkin(10) * t92 + t80;
t77 = pkin(5) * t100 - pkin(10) * t93 + t79;
t76 = t117 * t77 + t121 * t78;
t75 = -t117 * t78 + t121 * t77;
t1 = m(7) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t88 ^ 2) / 0.2e1 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,3) * t99 / 0.2e1) * t99 + (t88 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (-t88 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (t81 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,5) * t99 + Ifges(7,1) * t83 / 0.2e1) * t83 + (t95 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t102 / 0.2e1) * t102 + (-t81 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t83 + Ifges(7,6) * t99 + Ifges(7,2) * t82 / 0.2e1) * t82 + (-t95 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t102 + Ifges(5,6) * t114 + Ifges(5,2) * t101 / 0.2e1) * t101 + (t79 * mrSges(6,1) - t80 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t100 / 0.2e1) * t100 + (-t105 * mrSges(4,1) + t106 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * qJD(2)) * qJD(2) + ((-t104 * mrSges(4,1) + t106 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t124 + (t105 * mrSges(4,2) - t104 * mrSges(4,3) + (-pkin(7) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t120 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t120 ^ 2 + t124 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t130 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t124) * t124 + (-pkin(1) * mrSges(3,2) + (t130 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t120 + (Ifges(3,4) - Ifges(4,5)) * t124) * t120) * qJD(1)) * qJD(1);
T  = t1;
