% Calculate kinetic energy for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:33
% EndTime: 2019-03-08 23:15:34
% DurationCPUTime: 0.44s
% Computational Cost: add. (606->98), mult. (1309->156), div. (0->0), fcn. (961->12), ass. (0->45)
t131 = sin(qJ(2));
t125 = sin(pkin(6));
t141 = qJD(1) * t125;
t118 = qJD(2) * pkin(8) + t131 * t141;
t130 = sin(qJ(3));
t134 = cos(qJ(3));
t127 = cos(pkin(6));
t140 = qJD(1) * t127;
t112 = t134 * t118 + t130 * t140;
t108 = qJD(3) * pkin(9) + t112;
t135 = cos(qJ(2));
t137 = t135 * t141;
t113 = -t137 + (-pkin(3) * t134 - pkin(9) * t130 - pkin(2)) * qJD(2);
t129 = sin(qJ(4));
t133 = cos(qJ(4));
t102 = t133 * t108 + t129 * t113;
t139 = qJD(2) * t130;
t116 = t133 * qJD(3) - t129 * t139;
t100 = t116 * qJ(5) + t102;
t124 = sin(pkin(12));
t126 = cos(pkin(12));
t101 = -t129 * t108 + t133 * t113;
t117 = t129 * qJD(3) + t133 * t139;
t138 = t134 * qJD(2);
t122 = qJD(4) - t138;
t95 = t122 * pkin(4) - t117 * qJ(5) + t101;
t92 = t126 * t100 + t124 * t95;
t91 = -t124 * t100 + t126 * t95;
t111 = -t130 * t118 + t134 * t140;
t107 = -qJD(3) * pkin(3) - t111;
t103 = -t116 * pkin(4) + qJD(5) + t107;
t132 = cos(qJ(6));
t128 = sin(qJ(6));
t120 = qJD(6) + t122;
t119 = -qJD(2) * pkin(2) - t137;
t105 = t124 * t116 + t126 * t117;
t104 = t126 * t116 - t124 * t117;
t98 = t128 * t104 + t132 * t105;
t97 = t132 * t104 - t128 * t105;
t96 = -t104 * pkin(5) + t103;
t90 = t104 * pkin(10) + t92;
t89 = t122 * pkin(5) - t105 * pkin(10) + t91;
t88 = t128 * t89 + t132 * t90;
t87 = -t128 * t90 + t132 * t89;
t1 = m(6) * (t103 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t102 ^ 2 + t107 ^ 2) / 0.2e1 + (t96 * mrSges(7,2) - t87 * mrSges(7,3) + Ifges(7,1) * t98 / 0.2e1) * t98 + (t107 * mrSges(5,2) - t101 * mrSges(5,3) + Ifges(5,1) * t117 / 0.2e1) * t117 + (t103 * mrSges(6,2) - t91 * mrSges(6,3) + Ifges(6,1) * t105 / 0.2e1) * t105 + (t111 * mrSges(4,1) - t112 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t96 * mrSges(7,1) + t88 * mrSges(7,3) + Ifges(7,4) * t98 + Ifges(7,2) * t97 / 0.2e1) * t97 + (-t107 * mrSges(5,1) + t102 * mrSges(5,3) + Ifges(5,4) * t117 + Ifges(5,2) * t116 / 0.2e1) * t116 + (-t103 * mrSges(6,1) + t92 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,2) * t104 / 0.2e1) * t104 + (m(3) * (t127 ^ 2 + (t131 ^ 2 + t135 ^ 2) * t125 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t87 * mrSges(7,1) - t88 * mrSges(7,2) + Ifges(7,5) * t98 + Ifges(7,6) * t97 + Ifges(7,3) * t120 / 0.2e1) * t120 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t135 - mrSges(3,2) * t131) * t141 + (-t119 * mrSges(4,1) + t112 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t138 / 0.2e1) * t134 + (t119 * mrSges(4,2) - t111 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t134 + Ifges(4,1) * t130 / 0.2e1) * qJD(2)) * t130) * qJD(2) + (t101 * mrSges(5,1) + t91 * mrSges(6,1) - t102 * mrSges(5,2) - t92 * mrSges(6,2) + Ifges(5,5) * t117 + Ifges(6,5) * t105 + Ifges(5,6) * t116 + Ifges(6,6) * t104 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t122) * t122;
T  = t1;
