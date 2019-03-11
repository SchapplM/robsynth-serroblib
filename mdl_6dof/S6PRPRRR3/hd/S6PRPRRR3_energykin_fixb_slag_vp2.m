% Calculate kinetic energy for
% S6PRPRRR3
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:46
% EndTime: 2019-03-08 20:31:46
% DurationCPUTime: 0.40s
% Computational Cost: add. (533->93), mult. (1282->152), div. (0->0), fcn. (987->12), ass. (0->45)
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t124 = sin(pkin(12));
t126 = cos(pkin(12));
t130 = sin(qJ(4));
t134 = cos(qJ(4));
t116 = (t124 * t134 + t126 * t130) * qJD(2);
t131 = sin(qJ(2));
t125 = sin(pkin(6));
t139 = qJD(1) * t125;
t119 = qJD(2) * qJ(3) + t131 * t139;
t127 = cos(pkin(6));
t138 = qJD(1) * t127;
t121 = t126 * t138;
t140 = pkin(8) * qJD(2);
t110 = t121 + (-t119 - t140) * t124;
t113 = t126 * t119 + t124 * t138;
t111 = t126 * t140 + t113;
t99 = t134 * t110 - t111 * t130;
t97 = qJD(4) * pkin(4) - pkin(9) * t116 + t99;
t100 = t130 * t110 + t134 * t111;
t115 = (-t124 * t130 + t126 * t134) * qJD(2);
t98 = pkin(9) * t115 + t100;
t93 = t129 * t97 + t133 * t98;
t135 = cos(qJ(2));
t137 = -t135 * t139 + qJD(3);
t92 = -t129 * t98 + t133 * t97;
t104 = t115 * t133 - t116 * t129;
t114 = (-pkin(3) * t126 - pkin(2)) * qJD(2) + t137;
t106 = -pkin(4) * t115 + t114;
t132 = cos(qJ(6));
t128 = sin(qJ(6));
t123 = qJD(4) + qJD(5);
t118 = -qJD(2) * pkin(2) + t137;
t112 = -t119 * t124 + t121;
t105 = t115 * t129 + t116 * t133;
t103 = qJD(6) - t104;
t102 = t105 * t132 + t123 * t128;
t101 = -t105 * t128 + t123 * t132;
t94 = -pkin(5) * t104 - pkin(10) * t105 + t106;
t91 = pkin(10) * t123 + t93;
t90 = -pkin(5) * t123 - t92;
t89 = t128 * t94 + t132 * t91;
t88 = -t128 * t91 + t132 * t94;
t1 = m(7) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t114 ^ 2 + t99 ^ 2) / 0.2e1 + (t92 * mrSges(6,1) - t93 * mrSges(6,2) + Ifges(6,3) * t123 / 0.2e1) * t123 + (t114 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,1) * t116 / 0.2e1) * t116 + (t88 * mrSges(7,1) - t89 * mrSges(7,2) + Ifges(7,3) * t103 / 0.2e1) * t103 + (-t114 * mrSges(5,1) + t100 * mrSges(5,3) + Ifges(5,4) * t116 + Ifges(5,2) * t115 / 0.2e1) * t115 + (t106 * mrSges(6,2) - t92 * mrSges(6,3) + Ifges(6,5) * t123 + Ifges(6,1) * t105 / 0.2e1) * t105 + (t90 * mrSges(7,2) - t88 * mrSges(7,3) + Ifges(7,5) * t103 + Ifges(7,1) * t102 / 0.2e1) * t102 + (m(2) / 0.2e1 + m(3) * (t127 ^ 2 + (t131 ^ 2 + t135 ^ 2) * t125 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t106 * mrSges(6,1) + t93 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,6) * t123 + Ifges(6,2) * t104 / 0.2e1) * t104 + (-t90 * mrSges(7,1) + t89 * mrSges(7,3) + Ifges(7,4) * t102 + Ifges(7,6) * t103 + Ifges(7,2) * t101 / 0.2e1) * t101 + (t99 * mrSges(5,1) - t100 * mrSges(5,2) + Ifges(5,5) * t116 + Ifges(5,6) * t115 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t118 * (-mrSges(4,1) * t126 + mrSges(4,2) * t124) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t126 ^ 2 / 0.2e1 + (Ifges(4,4) * t126 + Ifges(4,1) * t124 / 0.2e1) * t124) * qJD(2) + (mrSges(3,1) * t135 - mrSges(3,2) * t131) * t139 + (-t112 * t124 + t113 * t126) * mrSges(4,3)) * qJD(2);
T  = t1;
