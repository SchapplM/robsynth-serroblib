% Calculate kinetic energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:45
% EndTime: 2019-03-09 22:38:46
% DurationCPUTime: 0.63s
% Computational Cost: add. (712->108), mult. (1445->156), div. (0->0), fcn. (1002->8), ass. (0->42)
t128 = -pkin(4) - pkin(5);
t127 = pkin(7) * mrSges(3,3);
t126 = cos(qJ(4));
t113 = sin(qJ(4));
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t115 = sin(qJ(2));
t124 = t115 * qJD(1);
t103 = qJD(2) * t114 + t117 * t124;
t118 = cos(qJ(2));
t125 = qJD(1) * t118;
t110 = qJD(3) - t125;
t101 = (-pkin(2) * t118 - pkin(8) * t115 - pkin(1)) * qJD(1);
t108 = pkin(7) * t125 + qJD(2) * pkin(8);
t95 = t117 * t101 - t108 * t114;
t89 = pkin(3) * t110 - pkin(9) * t103 + t95;
t102 = qJD(2) * t117 - t114 * t124;
t96 = t114 * t101 + t117 * t108;
t92 = pkin(9) * t102 + t96;
t84 = t113 * t89 + t126 * t92;
t109 = qJD(4) + t110;
t82 = t109 * qJ(5) + t84;
t107 = -qJD(2) * pkin(2) + pkin(7) * t124;
t83 = -t113 * t92 + t126 * t89;
t123 = qJD(5) - t83;
t122 = pkin(3) * t102 - t107;
t94 = t113 * t102 + t126 * t103;
t121 = qJ(5) * t94 + t122;
t116 = cos(qJ(6));
t112 = sin(qJ(6));
t106 = qJD(6) - t109;
t93 = -t126 * t102 + t103 * t113;
t87 = t112 * t93 + t116 * t94;
t86 = -t112 * t94 + t116 * t93;
t85 = pkin(4) * t93 - t121;
t81 = -t109 * pkin(4) + t123;
t80 = t128 * t93 + t121;
t79 = pkin(10) * t93 + t82;
t78 = -t94 * pkin(10) + t128 * t109 + t123;
t77 = t112 * t78 + t116 * t79;
t76 = -t112 * t79 + t116 * t78;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t107 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t122 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + (t80 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t87 / 0.2e1) * t87 + (t95 * mrSges(4,1) - t96 * mrSges(4,2) + Ifges(4,3) * t110 / 0.2e1) * t110 + (-t80 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t87 + Ifges(7,2) * t86 / 0.2e1) * t86 + (t107 * mrSges(4,2) - t95 * mrSges(4,3) + Ifges(4,5) * t110 + Ifges(4,1) * t103 / 0.2e1) * t103 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t106 / 0.2e1) * t106 + (-t107 * mrSges(4,1) + t96 * mrSges(4,3) + Ifges(4,4) * t103 + Ifges(4,6) * t110 + Ifges(4,2) * t102 / 0.2e1) * t102 + (-t122 * mrSges(5,2) + t81 * mrSges(6,2) - t83 * mrSges(5,3) - t85 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t94) * t94 + (-t122 * mrSges(5,1) + t85 * mrSges(6,1) - t82 * mrSges(6,2) - t84 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t93 + (-Ifges(5,4) + Ifges(6,5)) * t94) * t93 + (t83 * mrSges(5,1) - t81 * mrSges(6,1) - t84 * mrSges(5,2) + t82 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t109 + (Ifges(6,4) + Ifges(5,5)) * t94 + (-Ifges(5,6) + Ifges(6,6)) * t93) * t109 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t115 ^ 2 + t118 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t127) * t118) * t118 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t118 + (Ifges(3,1) / 0.2e1 + t127) * t115) * t115) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t118 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t115) * qJD(2)) * qJD(1);
T  = t1;
