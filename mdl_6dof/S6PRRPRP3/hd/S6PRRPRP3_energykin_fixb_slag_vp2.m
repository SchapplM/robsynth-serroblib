% Calculate kinetic energy for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:33:57
% EndTime: 2019-03-08 21:33:57
% DurationCPUTime: 0.40s
% Computational Cost: add. (430->94), mult. (951->139), div. (0->0), fcn. (661->10), ass. (0->38)
t127 = cos(qJ(5));
t116 = sin(qJ(5));
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t117 = sin(qJ(3));
t124 = qJD(2) * t117;
t106 = qJD(3) * t112 + t114 * t124;
t119 = cos(qJ(3));
t123 = qJD(2) * t119;
t120 = cos(qJ(2));
t113 = sin(pkin(6));
t126 = qJD(1) * t113;
t122 = t120 * t126;
t101 = -t122 + (-pkin(3) * t119 - qJ(4) * t117 - pkin(2)) * qJD(2);
t118 = sin(qJ(2));
t107 = qJD(2) * pkin(8) + t118 * t126;
t115 = cos(pkin(6));
t125 = qJD(1) * t115;
t100 = t119 * t107 + t117 * t125;
t98 = qJD(3) * qJ(4) + t100;
t89 = t114 * t101 - t112 * t98;
t86 = -pkin(4) * t123 - pkin(9) * t106 + t89;
t105 = qJD(3) * t114 - t112 * t124;
t90 = t112 * t101 + t114 * t98;
t88 = pkin(9) * t105 + t90;
t83 = t116 * t86 + t127 * t88;
t99 = -t117 * t107 + t119 * t125;
t82 = -t116 * t88 + t127 * t86;
t95 = -qJD(3) * pkin(3) + qJD(4) - t99;
t91 = -pkin(4) * t105 + t95;
t110 = qJD(5) - t123;
t108 = -qJD(2) * pkin(2) - t122;
t93 = t116 * t105 + t127 * t106;
t92 = -t127 * t105 + t106 * t116;
t84 = pkin(5) * t92 - qJ(6) * t93 + t91;
t81 = qJ(6) * t110 + t83;
t80 = -t110 * pkin(5) + qJD(6) - t82;
t1 = m(4) * (t100 ^ 2 + t108 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + (t95 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,1) * t106 / 0.2e1) * t106 + (t99 * mrSges(4,1) - t100 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t95 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,2) * t105 / 0.2e1) * t105 + (m(2) / 0.2e1 + m(3) * (t115 ^ 2 + (t118 ^ 2 + t120 ^ 2) * t113 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t120 - mrSges(3,2) * t118) * t126 + (t108 * mrSges(4,2) - t99 * mrSges(4,3) + Ifges(4,5) * qJD(3) + Ifges(4,1) * t124 / 0.2e1) * t117 + (-t108 * mrSges(4,1) - t89 * mrSges(5,1) + t90 * mrSges(5,2) + t100 * mrSges(4,3) - Ifges(5,5) * t106 + Ifges(4,6) * qJD(3) - Ifges(5,6) * t105 + (Ifges(4,4) * t117 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t119) * qJD(2)) * t119) * qJD(2) + (t91 * mrSges(6,2) + t80 * mrSges(7,2) - t82 * mrSges(6,3) - t84 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t93) * t93 + (t91 * mrSges(6,1) + t84 * mrSges(7,1) - t81 * mrSges(7,2) - t83 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t92 + (-Ifges(6,4) + Ifges(7,5)) * t93) * t92 + (t82 * mrSges(6,1) - t80 * mrSges(7,1) - t83 * mrSges(6,2) + t81 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t110 + (Ifges(7,4) + Ifges(6,5)) * t93 + (-Ifges(6,6) + Ifges(7,6)) * t92) * t110;
T  = t1;
