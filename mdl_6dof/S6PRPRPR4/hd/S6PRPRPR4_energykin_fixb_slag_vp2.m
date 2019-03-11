% Calculate kinetic energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:38:57
% EndTime: 2019-03-08 19:38:57
% DurationCPUTime: 0.39s
% Computational Cost: add. (501->93), mult. (1200->151), div. (0->0), fcn. (905->12), ass. (0->44)
t142 = cos(qJ(4));
t130 = cos(pkin(11));
t136 = cos(qJ(2));
t128 = sin(pkin(6));
t141 = qJD(1) * t128;
t138 = -t136 * t141 + qJD(3);
t115 = (-pkin(3) * t130 - pkin(2)) * qJD(2) + t138;
t127 = sin(pkin(11));
t133 = sin(qJ(4));
t139 = qJD(2) * t130;
t117 = qJD(2) * t127 * t133 - t139 * t142;
t118 = (t127 * t142 + t130 * t133) * qJD(2);
t105 = pkin(4) * t117 - qJ(5) * t118 + t115;
t126 = sin(pkin(12));
t129 = cos(pkin(12));
t134 = sin(qJ(2));
t121 = qJD(2) * qJ(3) + t134 * t141;
t131 = cos(pkin(6));
t140 = qJD(1) * t131;
t123 = t130 * t140;
t109 = t123 + (-pkin(8) * qJD(2) - t121) * t127;
t114 = t130 * t121 + t127 * t140;
t110 = pkin(8) * t139 + t114;
t100 = t133 * t109 + t110 * t142;
t98 = qJD(4) * qJ(5) + t100;
t94 = t126 * t105 + t129 * t98;
t93 = t129 * t105 - t126 * t98;
t99 = t109 * t142 - t133 * t110;
t97 = -qJD(4) * pkin(4) + qJD(5) - t99;
t135 = cos(qJ(6));
t132 = sin(qJ(6));
t120 = -qJD(2) * pkin(2) + t138;
t116 = qJD(6) + t117;
t113 = -t121 * t127 + t123;
t112 = qJD(4) * t126 + t118 * t129;
t111 = qJD(4) * t129 - t118 * t126;
t102 = t111 * t132 + t112 * t135;
t101 = t111 * t135 - t112 * t132;
t95 = -pkin(5) * t111 + t97;
t92 = pkin(9) * t111 + t94;
t91 = pkin(5) * t117 - pkin(9) * t112 + t93;
t90 = t132 * t91 + t135 * t92;
t89 = -t132 * t92 + t135 * t91;
t1 = m(4) * (t113 ^ 2 + t114 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t115 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t89 ^ 2 + t90 ^ 2 + t95 ^ 2) / 0.2e1 + (t115 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,1) * t118 / 0.2e1) * t118 + (t89 * mrSges(7,1) - t90 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t97 * mrSges(6,2) - t93 * mrSges(6,3) + Ifges(6,1) * t112 / 0.2e1) * t112 + (-t97 * mrSges(6,1) + t94 * mrSges(6,3) + Ifges(6,4) * t112 + Ifges(6,2) * t111 / 0.2e1) * t111 + (t95 * mrSges(7,2) - t89 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t102 / 0.2e1) * t102 + (t99 * mrSges(5,1) - t100 * mrSges(5,2) + Ifges(5,5) * t118 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (m(3) * (t131 ^ 2 + (t134 ^ 2 + t136 ^ 2) * t128 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t95 * mrSges(7,1) + t90 * mrSges(7,3) + Ifges(7,4) * t102 + Ifges(7,6) * t116 + Ifges(7,2) * t101 / 0.2e1) * t101 + (t115 * mrSges(5,1) + t93 * mrSges(6,1) - t94 * mrSges(6,2) - t100 * mrSges(5,3) - Ifges(5,4) * t118 + Ifges(6,5) * t112 - Ifges(5,6) * qJD(4) + Ifges(6,6) * t111 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t117) * t117 + (t120 * (-mrSges(4,1) * t130 + mrSges(4,2) * t127) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t130 ^ 2 / 0.2e1 + (Ifges(4,4) * t130 + Ifges(4,1) * t127 / 0.2e1) * t127) * qJD(2) + (mrSges(3,1) * t136 - mrSges(3,2) * t134) * t141 + (-t113 * t127 + t114 * t130) * mrSges(4,3)) * qJD(2);
T  = t1;
