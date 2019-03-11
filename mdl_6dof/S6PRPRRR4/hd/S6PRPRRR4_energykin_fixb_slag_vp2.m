% Calculate kinetic energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:50
% EndTime: 2019-03-08 20:35:51
% DurationCPUTime: 0.42s
% Computational Cost: add. (513->93), mult. (1200->153), div. (0->0), fcn. (905->12), ass. (0->45)
t130 = cos(pkin(12));
t139 = cos(qJ(2));
t129 = sin(pkin(6));
t144 = qJD(1) * t129;
t141 = -t139 * t144 + qJD(3);
t117 = (-pkin(3) * t130 - pkin(2)) * qJD(2) + t141;
t128 = sin(pkin(12));
t134 = sin(qJ(4));
t138 = cos(qJ(4));
t142 = qJD(2) * t130;
t119 = -t134 * t128 * qJD(2) + t138 * t142;
t120 = (t128 * t138 + t130 * t134) * qJD(2);
t106 = -pkin(4) * t119 - pkin(9) * t120 + t117;
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t135 = sin(qJ(2));
t123 = qJD(2) * qJ(3) + t135 * t144;
t131 = cos(pkin(6));
t143 = qJD(1) * t131;
t125 = t130 * t143;
t110 = t125 + (-pkin(8) * qJD(2) - t123) * t128;
t115 = t130 * t123 + t128 * t143;
t111 = pkin(8) * t142 + t115;
t101 = t134 * t110 + t138 * t111;
t99 = qJD(4) * pkin(9) + t101;
t95 = t133 * t106 + t137 * t99;
t94 = t137 * t106 - t133 * t99;
t100 = t110 * t138 - t134 * t111;
t118 = qJD(5) - t119;
t98 = -qJD(4) * pkin(4) - t100;
t136 = cos(qJ(6));
t132 = sin(qJ(6));
t122 = -qJD(2) * pkin(2) + t141;
t116 = qJD(6) + t118;
t114 = -t123 * t128 + t125;
t113 = qJD(4) * t133 + t120 * t137;
t112 = qJD(4) * t137 - t120 * t133;
t105 = t112 * t132 + t113 * t136;
t104 = t112 * t136 - t113 * t132;
t96 = -pkin(5) * t112 + t98;
t93 = pkin(10) * t112 + t95;
t92 = pkin(5) * t118 - pkin(10) * t113 + t94;
t91 = t132 * t92 + t136 * t93;
t90 = -t132 * t93 + t136 * t92;
t1 = m(5) * (t100 ^ 2 + t101 ^ 2 + t117 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t96 ^ 2) / 0.2e1 + (t117 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,1) * t120 / 0.2e1) * t120 + (t94 * mrSges(6,1) - t95 * mrSges(6,2) + Ifges(6,3) * t118 / 0.2e1) * t118 + (t90 * mrSges(7,1) - t91 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (-t117 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t120 + Ifges(5,2) * t119 / 0.2e1) * t119 + (t98 * mrSges(6,2) - t94 * mrSges(6,3) + Ifges(6,5) * t118 + Ifges(6,1) * t113 / 0.2e1) * t113 + (t96 * mrSges(7,2) - t90 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t105 / 0.2e1) * t105 + (m(3) * (t131 ^ 2 + (t135 ^ 2 + t139 ^ 2) * t129 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t98 * mrSges(6,1) + t95 * mrSges(6,3) + Ifges(6,4) * t113 + Ifges(6,6) * t118 + Ifges(6,2) * t112 / 0.2e1) * t112 + (-t96 * mrSges(7,1) + t91 * mrSges(7,3) + Ifges(7,4) * t105 + Ifges(7,6) * t116 + Ifges(7,2) * t104 / 0.2e1) * t104 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,5) * t120 + Ifges(5,6) * t119 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t122 * (-mrSges(4,1) * t130 + mrSges(4,2) * t128) + (Ifges(4,2) * t130 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t130 + Ifges(4,1) * t128 / 0.2e1) * t128) * qJD(2) + (mrSges(3,1) * t139 - mrSges(3,2) * t135) * t144 + (-t114 * t128 + t115 * t130) * mrSges(4,3)) * qJD(2);
T  = t1;
