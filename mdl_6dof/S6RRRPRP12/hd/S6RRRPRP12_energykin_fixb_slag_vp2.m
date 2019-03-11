% Calculate kinetic energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:32
% EndTime: 2019-03-09 17:51:33
% DurationCPUTime: 0.48s
% Computational Cost: add. (672->107), mult. (1476->149), div. (0->0), fcn. (1088->8), ass. (0->42)
t128 = pkin(3) + pkin(10);
t127 = cos(qJ(3));
t126 = cos(qJ(5));
t114 = sin(qJ(5));
t125 = cos(pkin(6)) * qJD(1);
t111 = qJD(2) + t125;
t115 = sin(qJ(3));
t116 = sin(qJ(2));
t112 = sin(pkin(6));
t124 = t112 * qJD(1);
t122 = t116 * t124;
t102 = t115 * t111 + t127 * t122;
t117 = cos(qJ(2));
t121 = t117 * t124;
t107 = qJD(3) - t121;
t123 = pkin(1) * t125;
t104 = pkin(8) * t121 + t116 * t123;
t97 = pkin(9) * t111 + t104;
t99 = (-pkin(2) * t117 - pkin(9) * t116 - pkin(1)) * t124;
t89 = -t115 * t97 + t127 * t99;
t120 = qJD(4) - t89;
t82 = t102 * pkin(4) - t128 * t107 + t120;
t101 = -t127 * t111 + t115 * t122;
t103 = -pkin(8) * t122 + t117 * t123;
t96 = -pkin(2) * t111 - t103;
t119 = -qJ(4) * t102 + t96;
t84 = t128 * t101 + t119;
t79 = t114 * t82 + t126 * t84;
t90 = t115 * t99 + t127 * t97;
t88 = -t107 * qJ(4) - t90;
t85 = -pkin(4) * t101 - t88;
t78 = -t114 * t84 + t126 * t82;
t118 = qJD(1) ^ 2;
t100 = qJD(5) + t102;
t92 = t114 * t101 + t126 * t107;
t91 = -t126 * t101 + t107 * t114;
t87 = -t107 * pkin(3) + t120;
t86 = pkin(3) * t101 + t119;
t80 = pkin(5) * t91 - qJ(6) * t92 + t85;
t77 = qJ(6) * t100 + t79;
t76 = -t100 * pkin(5) + qJD(6) - t78;
t1 = m(4) * (t89 ^ 2 + t90 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t85 ^ 2) / 0.2e1 + t118 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t112 ^ 2 * t118 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t117 / 0.2e1) * t117 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t117 + Ifges(3,1) * t116 / 0.2e1) * t116) * t124 + (-t103 * t116 + t104 * t117) * mrSges(3,3)) * t124 + (t103 * mrSges(3,1) - t104 * mrSges(3,2) + Ifges(3,3) * t111 / 0.2e1 + (Ifges(3,5) * t116 + Ifges(3,6) * t117) * t124) * t111 + (t85 * mrSges(6,2) + t76 * mrSges(7,2) - t78 * mrSges(6,3) - t80 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t92) * t92 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + t87 * mrSges(5,2) - t88 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t107) * t107 + (t85 * mrSges(6,1) + t80 * mrSges(7,1) - t77 * mrSges(7,2) - t79 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t91 + (-Ifges(6,4) + Ifges(7,5)) * t92) * t91 + (t87 * mrSges(5,1) + t96 * mrSges(4,2) - t89 * mrSges(4,3) - t86 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t102 + (-Ifges(5,4) + Ifges(4,5)) * t107) * t102 + (t96 * mrSges(4,1) + t88 * mrSges(5,1) - t86 * mrSges(5,2) - t90 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101 + (Ifges(5,5) - Ifges(4,6)) * t107 + (-Ifges(4,4) - Ifges(5,6)) * t102) * t101 + (t78 * mrSges(6,1) - t76 * mrSges(7,1) - t79 * mrSges(6,2) + t77 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t100 + (Ifges(7,4) + Ifges(6,5)) * t92 + (-Ifges(6,6) + Ifges(7,6)) * t91) * t100;
T  = t1;
