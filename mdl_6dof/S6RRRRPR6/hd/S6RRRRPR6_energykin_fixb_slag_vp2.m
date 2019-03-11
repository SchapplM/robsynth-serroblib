% Calculate kinetic energy for
% S6RRRRPR6
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:34
% EndTime: 2019-03-09 22:18:34
% DurationCPUTime: 0.64s
% Computational Cost: add. (1166->111), mult. (2493->171), div. (0->0), fcn. (1874->10), ass. (0->45)
t131 = pkin(7) * mrSges(3,3);
t117 = sin(pkin(11));
t118 = cos(pkin(11));
t121 = sin(qJ(3));
t125 = cos(qJ(3));
t122 = sin(qJ(2));
t130 = qJD(1) * t122;
t108 = qJD(2) * t125 - t121 * t130;
t109 = qJD(2) * t121 + t125 * t130;
t120 = sin(qJ(4));
t124 = cos(qJ(4));
t101 = t108 * t120 + t109 * t124;
t126 = cos(qJ(2));
t129 = qJD(1) * t126;
t115 = qJD(3) - t129;
t114 = qJD(4) + t115;
t107 = (-pkin(2) * t126 - pkin(8) * t122 - pkin(1)) * qJD(1);
t113 = pkin(7) * t129 + qJD(2) * pkin(8);
t102 = t125 * t107 - t113 * t121;
t97 = pkin(3) * t115 - pkin(9) * t109 + t102;
t103 = t121 * t107 + t125 * t113;
t99 = pkin(9) * t108 + t103;
t90 = -t120 * t99 + t124 * t97;
t86 = pkin(4) * t114 - qJ(5) * t101 + t90;
t100 = t108 * t124 - t109 * t120;
t91 = t120 * t97 + t124 * t99;
t88 = qJ(5) * t100 + t91;
t81 = t117 * t86 + t118 * t88;
t80 = -t117 * t88 + t118 * t86;
t112 = -qJD(2) * pkin(2) + pkin(7) * t130;
t104 = -pkin(3) * t108 + t112;
t94 = -pkin(4) * t100 + qJD(5) + t104;
t123 = cos(qJ(6));
t119 = sin(qJ(6));
t111 = qJD(6) + t114;
t93 = t100 * t117 + t101 * t118;
t92 = t100 * t118 - t101 * t117;
t89 = -pkin(5) * t92 + t94;
t85 = t119 * t92 + t123 * t93;
t84 = -t119 * t93 + t123 * t92;
t79 = pkin(10) * t92 + t81;
t78 = pkin(5) * t114 - pkin(10) * t93 + t80;
t77 = t119 * t78 + t123 * t79;
t76 = -t119 * t79 + t123 * t78;
t1 = m(4) * (t102 ^ 2 + t103 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t76 ^ 2 + t77 ^ 2 + t89 ^ 2) / 0.2e1 + (t94 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t89 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t85 / 0.2e1) * t85 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,3) * t115 / 0.2e1) * t115 + (t104 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,1) * t101 / 0.2e1) * t101 + (-t94 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t89 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t112 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,5) * t115 + Ifges(4,1) * t109 / 0.2e1) * t109 + (-t104 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t101 + Ifges(5,2) * t100 / 0.2e1) * t100 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t85 + Ifges(7,6) * t84 + Ifges(7,3) * t111 / 0.2e1) * t111 + (-t112 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t109 + Ifges(4,6) * t115 + Ifges(4,2) * t108 / 0.2e1) * t108 + (t90 * mrSges(5,1) + t80 * mrSges(6,1) - t91 * mrSges(5,2) - t81 * mrSges(6,2) + Ifges(5,5) * t101 + Ifges(6,5) * t93 + Ifges(5,6) * t100 + Ifges(6,6) * t92 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t114) * t114 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t122 ^ 2 + t126 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t131) * t126) * t126 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t126 + (Ifges(3,1) / 0.2e1 + t131) * t122) * t122) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t126 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t122) * qJD(2)) * qJD(1);
T  = t1;
