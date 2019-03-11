% Calculate kinetic energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:26
% EndTime: 2019-03-09 22:26:27
% DurationCPUTime: 0.67s
% Computational Cost: add. (1488->114), mult. (3378->181), div. (0->0), fcn. (2746->12), ass. (0->52)
t131 = sin(pkin(12));
t133 = cos(pkin(12));
t147 = cos(pkin(6)) * qJD(1);
t130 = qJD(2) + t147;
t137 = sin(qJ(3));
t141 = cos(qJ(3));
t138 = sin(qJ(2));
t132 = sin(pkin(6));
t148 = qJD(1) * t132;
t145 = t138 * t148;
t121 = t130 * t141 - t137 * t145;
t122 = t130 * t137 + t141 * t145;
t136 = sin(qJ(4));
t140 = cos(qJ(4));
t113 = t121 * t136 + t122 * t140;
t142 = cos(qJ(2));
t144 = t142 * t148;
t126 = qJD(3) - t144;
t125 = qJD(4) + t126;
t146 = pkin(1) * t147;
t124 = pkin(8) * t144 + t138 * t146;
t119 = pkin(9) * t130 + t124;
t120 = (-pkin(2) * t142 - pkin(9) * t138 - pkin(1)) * t148;
t110 = -t119 * t137 + t141 * t120;
t107 = pkin(3) * t126 - pkin(10) * t122 + t110;
t111 = t141 * t119 + t137 * t120;
t109 = pkin(10) * t121 + t111;
t97 = t140 * t107 - t109 * t136;
t94 = pkin(4) * t125 - qJ(5) * t113 + t97;
t112 = t121 * t140 - t122 * t136;
t98 = t136 * t107 + t140 * t109;
t96 = qJ(5) * t112 + t98;
t91 = t131 * t94 + t133 * t96;
t90 = -t131 * t96 + t133 * t94;
t102 = t112 * t133 - t113 * t131;
t123 = -pkin(8) * t145 + t142 * t146;
t118 = -pkin(2) * t130 - t123;
t114 = -pkin(3) * t121 + t118;
t104 = -pkin(4) * t112 + qJD(5) + t114;
t143 = qJD(1) ^ 2;
t139 = cos(qJ(6));
t135 = sin(qJ(6));
t103 = t112 * t131 + t113 * t133;
t101 = qJD(6) - t102;
t100 = t103 * t139 + t125 * t135;
t99 = -t103 * t135 + t125 * t139;
t92 = -pkin(5) * t102 - pkin(11) * t103 + t104;
t89 = pkin(11) * t125 + t91;
t88 = -pkin(5) * t125 - t90;
t87 = t135 * t92 + t139 * t89;
t86 = -t135 * t89 + t139 * t92;
t1 = m(6) * (t104 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t132 ^ 2 * t143 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + t143 * Ifges(2,3) / 0.2e1 + m(5) * (t114 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t118 ^ 2) / 0.2e1 + (-t88 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t110 * mrSges(4,1) - t111 * mrSges(4,2) + Ifges(4,3) * t126 / 0.2e1) * t126 + (t114 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,1) * t113 / 0.2e1) * t113 + (t104 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,1) * t103 / 0.2e1) * t103 + (t118 * mrSges(4,2) - t110 * mrSges(4,3) + Ifges(4,5) * t126 + Ifges(4,1) * t122 / 0.2e1) * t122 + (-t114 * mrSges(5,1) + t98 * mrSges(5,3) + Ifges(5,4) * t113 + Ifges(5,2) * t112 / 0.2e1) * t112 + (-t104 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t103 + Ifges(6,2) * t102 / 0.2e1) * t102 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t101 / 0.2e1) * t101 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t142 / 0.2e1) * t142 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t142 + Ifges(3,1) * t138 / 0.2e1) * t138) * t148 + (-t123 * t138 + t124 * t142) * mrSges(3,3)) * t148 + (t123 * mrSges(3,1) - t124 * mrSges(3,2) + Ifges(3,3) * t130 / 0.2e1 + (Ifges(3,5) * t138 + Ifges(3,6) * t142) * t148) * t130 + (-t118 * mrSges(4,1) + t111 * mrSges(4,3) + Ifges(4,4) * t122 + Ifges(4,6) * t126 + Ifges(4,2) * t121 / 0.2e1) * t121 + (t88 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t101 + Ifges(7,1) * t100 / 0.2e1) * t100 + (t97 * mrSges(5,1) + t90 * mrSges(6,1) - t98 * mrSges(5,2) - t91 * mrSges(6,2) + Ifges(5,5) * t113 + Ifges(6,5) * t103 + Ifges(5,6) * t112 + Ifges(6,6) * t102 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t125) * t125;
T  = t1;
