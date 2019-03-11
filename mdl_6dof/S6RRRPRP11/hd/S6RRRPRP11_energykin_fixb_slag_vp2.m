% Calculate kinetic energy for
% S6RRRPRP11
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:19
% EndTime: 2019-03-09 17:41:19
% DurationCPUTime: 0.46s
% Computational Cost: add. (674->107), mult. (1484->149), div. (0->0), fcn. (1096->8), ass. (0->42)
t128 = pkin(3) + pkin(10);
t127 = cos(qJ(3));
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t125 = cos(pkin(6)) * qJD(1);
t111 = qJD(2) + t125;
t115 = sin(qJ(3));
t116 = sin(qJ(2));
t112 = sin(pkin(6));
t126 = qJD(1) * t112;
t123 = t116 * t126;
t102 = t115 * t111 + t127 * t123;
t118 = cos(qJ(2));
t122 = t118 * t126;
t107 = qJD(3) - t122;
t124 = pkin(1) * t125;
t104 = pkin(8) * t122 + t116 * t124;
t98 = pkin(9) * t111 + t104;
t99 = (-pkin(2) * t118 - pkin(9) * t116 - pkin(1)) * t126;
t90 = -t115 * t98 + t127 * t99;
t121 = qJD(4) - t90;
t83 = t102 * pkin(4) - t128 * t107 + t121;
t101 = -t127 * t111 + t115 * t123;
t103 = -pkin(8) * t123 + t118 * t124;
t97 = -pkin(2) * t111 - t103;
t120 = -qJ(4) * t102 + t97;
t85 = t128 * t101 + t120;
t79 = t114 * t83 + t117 * t85;
t91 = t115 * t99 + t127 * t98;
t89 = -t107 * qJ(4) - t91;
t78 = -t114 * t85 + t117 * t83;
t86 = -pkin(4) * t101 - t89;
t119 = qJD(1) ^ 2;
t100 = qJD(5) + t102;
t93 = t101 * t114 + t107 * t117;
t92 = t101 * t117 - t107 * t114;
t88 = -t107 * pkin(3) + t121;
t87 = pkin(3) * t101 + t120;
t80 = -pkin(5) * t92 + qJD(6) + t86;
t77 = qJ(6) * t92 + t79;
t76 = pkin(5) * t100 - qJ(6) * t93 + t78;
t1 = m(3) * (pkin(1) ^ 2 * t112 ^ 2 * t119 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t119 * Ifges(2,3) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t118 / 0.2e1) * t118 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t118 + Ifges(3,1) * t116 / 0.2e1) * t116) * t126 + (-t103 * t116 + t104 * t118) * mrSges(3,3)) * t126 + (t103 * mrSges(3,1) - t104 * mrSges(3,2) + Ifges(3,3) * t111 / 0.2e1 + (Ifges(3,5) * t116 + Ifges(3,6) * t118) * t126) * t111 + (t86 * mrSges(6,2) + t80 * mrSges(7,2) - t78 * mrSges(6,3) - t76 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t93) * t93 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + t88 * mrSges(5,2) - t89 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t107) * t107 + (-t86 * mrSges(6,1) - t80 * mrSges(7,1) + t79 * mrSges(6,3) + t77 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t92 + (Ifges(6,4) + Ifges(7,4)) * t93) * t92 + (t88 * mrSges(5,1) + t97 * mrSges(4,2) - t90 * mrSges(4,3) - t87 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t102 + (-Ifges(5,4) + Ifges(4,5)) * t107) * t102 + (t97 * mrSges(4,1) + t89 * mrSges(5,1) - t87 * mrSges(5,2) - t91 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101 + (Ifges(5,5) - Ifges(4,6)) * t107 + (-Ifges(4,4) - Ifges(5,6)) * t102) * t101 + (t78 * mrSges(6,1) + t76 * mrSges(7,1) - t79 * mrSges(6,2) - t77 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t100 + (Ifges(6,5) + Ifges(7,5)) * t93 + (Ifges(6,6) + Ifges(7,6)) * t92) * t100;
T  = t1;
