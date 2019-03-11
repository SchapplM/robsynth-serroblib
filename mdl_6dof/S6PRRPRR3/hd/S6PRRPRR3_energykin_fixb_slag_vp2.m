% Calculate kinetic energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:18
% EndTime: 2019-03-08 22:02:18
% DurationCPUTime: 0.53s
% Computational Cost: add. (710->102), mult. (1869->168), div. (0->0), fcn. (1509->14), ass. (0->51)
t149 = cos(qJ(2));
t138 = sin(pkin(6));
t156 = qJD(1) * t138;
t130 = qJD(2) * pkin(2) + t149 * t156;
t137 = sin(pkin(7));
t140 = cos(pkin(7));
t141 = cos(pkin(6));
t155 = qJD(1) * t141;
t159 = t130 * t140 + t137 * t155;
t136 = sin(pkin(13));
t139 = cos(pkin(13));
t144 = sin(qJ(3));
t148 = cos(qJ(3));
t154 = qJD(2) * t137;
t124 = (t136 * t144 - t139 * t148) * t154;
t145 = sin(qJ(2));
t129 = pkin(9) * t154 + t145 * t156;
t134 = qJD(2) * t140 + qJD(3);
t153 = qJ(4) * t154;
t157 = t159 * t148;
t112 = pkin(3) * t134 + (-t129 - t153) * t144 + t157;
t117 = t148 * t129 + t144 * t159;
t115 = t148 * t153 + t117;
t106 = t112 * t136 + t115 * t139;
t104 = pkin(10) * t134 + t106;
t133 = t140 * t155;
t121 = qJD(4) + t133 + (-pkin(3) * qJD(2) * t148 - t130) * t137;
t125 = (t136 * t148 + t139 * t144) * t154;
t108 = pkin(4) * t124 - pkin(10) * t125 + t121;
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t100 = t104 * t147 + t108 * t143;
t105 = t112 * t139 - t115 * t136;
t99 = -t104 * t143 + t108 * t147;
t119 = -t125 * t143 + t134 * t147;
t103 = -pkin(4) * t134 - t105;
t146 = cos(qJ(6));
t142 = sin(qJ(6));
t123 = qJD(5) + t124;
t122 = -t130 * t137 + t133;
t120 = t125 * t147 + t134 * t143;
t118 = qJD(6) - t119;
t116 = -t129 * t144 + t157;
t111 = t120 * t146 + t123 * t142;
t110 = -t120 * t142 + t123 * t146;
t101 = -pkin(5) * t119 - pkin(11) * t120 + t103;
t98 = pkin(11) * t123 + t100;
t97 = -pkin(5) * t123 - t99;
t96 = t101 * t142 + t146 * t98;
t95 = t101 * t146 - t142 * t98;
t1 = m(5) * (t105 ^ 2 + t106 ^ 2 + t121 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + (t121 * mrSges(5,2) - t105 * mrSges(5,3) + Ifges(5,1) * t125 / 0.2e1) * t125 + (t99 * mrSges(6,1) - t100 * mrSges(6,2) + Ifges(6,3) * t123 / 0.2e1) * t123 + (t95 * mrSges(7,1) - t96 * mrSges(7,2) + Ifges(7,3) * t118 / 0.2e1) * t118 - (-t121 * mrSges(5,1) + t106 * mrSges(5,3) + Ifges(5,4) * t125 - Ifges(5,2) * t124 / 0.2e1) * t124 + (t103 * mrSges(6,2) - t99 * mrSges(6,3) + Ifges(6,5) * t123 + Ifges(6,1) * t120 / 0.2e1) * t120 + (t97 * mrSges(7,2) - t95 * mrSges(7,3) + Ifges(7,5) * t118 + Ifges(7,1) * t111 / 0.2e1) * t111 + (m(3) * (t141 ^ 2 + (t145 ^ 2 + t149 ^ 2) * t138 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t103 * mrSges(6,1) + t100 * mrSges(6,3) + Ifges(6,4) * t120 + Ifges(6,6) * t123 + Ifges(6,2) * t119 / 0.2e1) * t119 + (-t97 * mrSges(7,1) + t96 * mrSges(7,3) + Ifges(7,4) * t111 + Ifges(7,6) * t118 + Ifges(7,2) * t110 / 0.2e1) * t110 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t149 - mrSges(3,2) * t145) * t156 + (t122 * (-mrSges(4,1) * t148 + mrSges(4,2) * t144) + (Ifges(4,2) * t148 ^ 2 / 0.2e1 + (Ifges(4,4) * t148 + Ifges(4,1) * t144 / 0.2e1) * t144) * t154 + (-t116 * t144 + t117 * t148) * mrSges(4,3)) * t137) * qJD(2) + (t116 * mrSges(4,1) + t105 * mrSges(5,1) - t117 * mrSges(4,2) - t106 * mrSges(5,2) + Ifges(5,5) * t125 - Ifges(5,6) * t124 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t134 + (Ifges(4,5) * t144 + Ifges(4,6) * t148) * t154) * t134;
T  = t1;
