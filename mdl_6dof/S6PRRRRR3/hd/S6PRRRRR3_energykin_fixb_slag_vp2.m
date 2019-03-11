% Calculate kinetic energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:57
% EndTime: 2019-03-09 00:47:57
% DurationCPUTime: 0.46s
% Computational Cost: add. (618->98), mult. (1309->158), div. (0->0), fcn. (961->12), ass. (0->46)
t132 = sin(qJ(2));
t126 = sin(pkin(6));
t143 = qJD(1) * t126;
t119 = qJD(2) * pkin(8) + t132 * t143;
t131 = sin(qJ(3));
t136 = cos(qJ(3));
t127 = cos(pkin(6));
t142 = qJD(1) * t127;
t113 = t136 * t119 + t131 * t142;
t109 = qJD(3) * pkin(9) + t113;
t137 = cos(qJ(2));
t139 = t137 * t143;
t114 = -t139 + (-pkin(3) * t136 - pkin(9) * t131 - pkin(2)) * qJD(2);
t130 = sin(qJ(4));
t135 = cos(qJ(4));
t103 = t135 * t109 + t130 * t114;
t141 = qJD(2) * t131;
t117 = qJD(3) * t135 - t130 * t141;
t101 = pkin(10) * t117 + t103;
t129 = sin(qJ(5));
t134 = cos(qJ(5));
t102 = -t109 * t130 + t135 * t114;
t118 = qJD(3) * t130 + t135 * t141;
t140 = qJD(2) * t136;
t124 = qJD(4) - t140;
t96 = pkin(4) * t124 - pkin(10) * t118 + t102;
t93 = t134 * t101 + t129 * t96;
t92 = -t101 * t129 + t134 * t96;
t112 = -t131 * t119 + t136 * t142;
t122 = qJD(5) + t124;
t108 = -qJD(3) * pkin(3) - t112;
t104 = -pkin(4) * t117 + t108;
t133 = cos(qJ(6));
t128 = sin(qJ(6));
t121 = qJD(6) + t122;
t120 = -qJD(2) * pkin(2) - t139;
t106 = t117 * t129 + t118 * t134;
t105 = t117 * t134 - t118 * t129;
t100 = t105 * t128 + t106 * t133;
t99 = t105 * t133 - t106 * t128;
t97 = -pkin(5) * t105 + t104;
t91 = pkin(11) * t105 + t93;
t90 = pkin(5) * t122 - pkin(11) * t106 + t92;
t89 = t128 * t90 + t133 * t91;
t88 = -t128 * t91 + t133 * t90;
t1 = m(7) * (t88 ^ 2 + t89 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t108 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t120 ^ 2) / 0.2e1 + (-t97 * mrSges(7,1) + t89 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t102 * mrSges(5,1) - t103 * mrSges(5,2) + Ifges(5,3) * t124 / 0.2e1) * t124 + (t92 * mrSges(6,1) - t93 * mrSges(6,2) + Ifges(6,3) * t122 / 0.2e1) * t122 + (t112 * mrSges(4,1) - t113 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t88 * mrSges(7,1) - t89 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t121 / 0.2e1) * t121 + (t108 * mrSges(5,2) - t102 * mrSges(5,3) + Ifges(5,5) * t124 + Ifges(5,1) * t118 / 0.2e1) * t118 + (t104 * mrSges(6,2) - t92 * mrSges(6,3) + Ifges(6,5) * t122 + Ifges(6,1) * t106 / 0.2e1) * t106 + (m(3) * (t127 ^ 2 + (t132 ^ 2 + t137 ^ 2) * t126 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t108 * mrSges(5,1) + t103 * mrSges(5,3) + Ifges(5,4) * t118 + Ifges(5,6) * t124 + Ifges(5,2) * t117 / 0.2e1) * t117 + (-t104 * mrSges(6,1) + t93 * mrSges(6,3) + Ifges(6,4) * t106 + Ifges(6,6) * t122 + Ifges(6,2) * t105 / 0.2e1) * t105 + (t97 * mrSges(7,2) - t88 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t121 + Ifges(7,1) * t100 / 0.2e1) * t100 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t137 - mrSges(3,2) * t132) * t143 + (-t120 * mrSges(4,1) + t113 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t140 / 0.2e1) * t136 + (t120 * mrSges(4,2) - t112 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t136 + Ifges(4,1) * t131 / 0.2e1) * qJD(2)) * t131) * qJD(2);
T  = t1;
