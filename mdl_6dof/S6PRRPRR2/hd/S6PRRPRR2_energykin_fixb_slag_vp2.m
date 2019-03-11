% Calculate kinetic energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:55:58
% EndTime: 2019-03-08 21:55:59
% DurationCPUTime: 0.43s
% Computational Cost: add. (538->99), mult. (1231->160), div. (0->0), fcn. (907->12), ass. (0->45)
t135 = cos(qJ(3));
t136 = cos(qJ(2));
t126 = sin(pkin(6));
t141 = qJD(1) * t126;
t138 = t136 * t141;
t114 = -t138 + qJD(4) + (-pkin(3) * t135 - pkin(2)) * qJD(2);
t125 = sin(pkin(12));
t127 = cos(pkin(12));
t131 = sin(qJ(3));
t139 = qJD(2) * t135;
t116 = -t125 * t131 * qJD(2) + t127 * t139;
t117 = (t125 * t135 + t127 * t131) * qJD(2);
t103 = -pkin(4) * t116 - pkin(9) * t117 + t114;
t130 = sin(qJ(5));
t134 = cos(qJ(5));
t132 = sin(qJ(2));
t119 = qJD(2) * pkin(8) + t132 * t141;
t128 = cos(pkin(6));
t140 = qJD(1) * t128;
t123 = t135 * t140;
t107 = qJD(3) * pkin(3) + t123 + (-qJ(4) * qJD(2) - t119) * t131;
t112 = t135 * t119 + t131 * t140;
t108 = qJ(4) * t139 + t112;
t98 = t125 * t107 + t127 * t108;
t96 = qJD(3) * pkin(9) + t98;
t92 = t130 * t103 + t134 * t96;
t91 = t134 * t103 - t130 * t96;
t97 = t107 * t127 - t125 * t108;
t115 = qJD(5) - t116;
t95 = -qJD(3) * pkin(4) - t97;
t133 = cos(qJ(6));
t129 = sin(qJ(6));
t120 = -qJD(2) * pkin(2) - t138;
t113 = qJD(6) + t115;
t111 = -t119 * t131 + t123;
t110 = qJD(3) * t130 + t117 * t134;
t109 = qJD(3) * t134 - t117 * t130;
t100 = t109 * t129 + t110 * t133;
t99 = t109 * t133 - t110 * t129;
t93 = -pkin(5) * t109 + t95;
t90 = pkin(10) * t109 + t92;
t89 = pkin(5) * t115 - pkin(10) * t110 + t91;
t88 = t129 * t89 + t133 * t90;
t87 = -t129 * t90 + t133 * t89;
t1 = m(5) * (t114 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t120 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t93 ^ 2) / 0.2e1 + (-t93 * mrSges(7,1) + t88 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t114 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,1) * t117 / 0.2e1) * t117 + (t91 * mrSges(6,1) - t92 * mrSges(6,2) + Ifges(6,3) * t115 / 0.2e1) * t115 + (-t114 * mrSges(5,1) + t98 * mrSges(5,3) + Ifges(5,4) * t117 + Ifges(5,2) * t116 / 0.2e1) * t116 + (t87 * mrSges(7,1) - t88 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t113 / 0.2e1) * t113 + (t95 * mrSges(6,2) - t91 * mrSges(6,3) + Ifges(6,5) * t115 + Ifges(6,1) * t110 / 0.2e1) * t110 + (m(3) * (t128 ^ 2 + (t132 ^ 2 + t136 ^ 2) * t126 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t95 * mrSges(6,1) + t92 * mrSges(6,3) + Ifges(6,4) * t110 + Ifges(6,6) * t115 + Ifges(6,2) * t109 / 0.2e1) * t109 + (t93 * mrSges(7,2) - t87 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t113 + Ifges(7,1) * t100 / 0.2e1) * t100 + (t120 * (-mrSges(4,1) * t135 + mrSges(4,2) * t131) + (Ifges(4,2) * t135 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t135 + Ifges(4,1) * t131 / 0.2e1) * t131) * qJD(2) + (mrSges(3,1) * t136 - mrSges(3,2) * t132) * t141 + (-t111 * t131 + t112 * t135) * mrSges(4,3)) * qJD(2) + (t111 * mrSges(4,1) + t97 * mrSges(5,1) - t112 * mrSges(4,2) - t98 * mrSges(5,2) + Ifges(5,5) * t117 + Ifges(5,6) * t116 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t131 + Ifges(4,6) * t135) * qJD(2)) * qJD(3);
T  = t1;
