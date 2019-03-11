% Calculate kinetic energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:13
% EndTime: 2019-03-10 04:48:14
% DurationCPUTime: 0.86s
% Computational Cost: add. (2214->122), mult. (5644->199), div. (0->0), fcn. (4776->14), ass. (0->60)
t156 = sin(qJ(2));
t161 = cos(qJ(2));
t149 = sin(pkin(6));
t169 = qJD(1) * t149;
t165 = t161 * t169;
t168 = cos(pkin(6)) * qJD(1);
t167 = pkin(1) * t168;
t142 = pkin(9) * t165 + t156 * t167;
t147 = qJD(2) + t168;
t148 = sin(pkin(7));
t150 = cos(pkin(7));
t163 = t147 * t148 + t150 * t165;
t131 = pkin(10) * t163 + t142;
t146 = t161 * t167;
t166 = t156 * t169;
t134 = pkin(2) * t147 + t146 + (-pkin(10) * t150 - pkin(9)) * t166;
t139 = (-pkin(10) * t148 * t156 - pkin(2) * t161 - pkin(1)) * t169;
t155 = sin(qJ(3));
t160 = cos(qJ(3));
t123 = -t155 * t131 + (t134 * t150 + t139 * t148) * t160;
t135 = -t155 * t166 + t160 * t163;
t171 = t148 * t155;
t170 = t150 * t155;
t125 = -t134 * t148 + t150 * t139;
t136 = t147 * t171 + (t156 * t160 + t161 * t170) * t169;
t116 = -pkin(3) * t135 - pkin(11) * t136 + t125;
t124 = t160 * t131 + t134 * t170 + t139 * t171;
t140 = t147 * t150 - t148 * t165 + qJD(3);
t119 = pkin(11) * t140 + t124;
t154 = sin(qJ(4));
t159 = cos(qJ(4));
t109 = t159 * t116 - t119 * t154;
t127 = t136 * t159 + t140 * t154;
t133 = qJD(4) - t135;
t106 = pkin(4) * t133 - pkin(12) * t127 + t109;
t110 = t154 * t116 + t159 * t119;
t126 = -t136 * t154 + t140 * t159;
t108 = pkin(12) * t126 + t110;
t153 = sin(qJ(5));
t158 = cos(qJ(5));
t103 = t153 * t106 + t158 * t108;
t102 = t106 * t158 - t108 * t153;
t121 = t126 * t158 - t127 * t153;
t118 = -pkin(3) * t140 - t123;
t111 = -pkin(4) * t126 + t118;
t162 = qJD(1) ^ 2;
t157 = cos(qJ(6));
t152 = sin(qJ(6));
t141 = -pkin(9) * t166 + t146;
t132 = qJD(5) + t133;
t122 = t126 * t153 + t127 * t158;
t120 = qJD(6) - t121;
t113 = t122 * t157 + t132 * t152;
t112 = -t122 * t152 + t132 * t157;
t104 = -pkin(5) * t121 - pkin(13) * t122 + t111;
t101 = pkin(13) * t132 + t103;
t100 = -pkin(5) * t132 - t102;
t99 = t101 * t157 + t104 * t152;
t98 = -t101 * t152 + t104 * t157;
t1 = m(5) * (t109 ^ 2 + t110 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + t162 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t149 ^ 2 * t162 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + (t123 * mrSges(4,1) - t124 * mrSges(4,2) + Ifges(4,3) * t140 / 0.2e1) * t140 + (t109 * mrSges(5,1) - t110 * mrSges(5,2) + Ifges(5,3) * t133 / 0.2e1) * t133 + (t102 * mrSges(6,1) - t103 * mrSges(6,2) + Ifges(6,3) * t132 / 0.2e1) * t132 + (t98 * mrSges(7,1) - t99 * mrSges(7,2) + Ifges(7,3) * t120 / 0.2e1) * t120 + (t125 * mrSges(4,2) - t123 * mrSges(4,3) + Ifges(4,5) * t140 + Ifges(4,1) * t136 / 0.2e1) * t136 + (t118 * mrSges(5,2) - t109 * mrSges(5,3) + Ifges(5,5) * t133 + Ifges(5,1) * t127 / 0.2e1) * t127 + (t111 * mrSges(6,2) - t102 * mrSges(6,3) + Ifges(6,5) * t132 + Ifges(6,1) * t122 / 0.2e1) * t122 + (t100 * mrSges(7,2) - t98 * mrSges(7,3) + Ifges(7,5) * t120 + Ifges(7,1) * t113 / 0.2e1) * t113 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t161 / 0.2e1) * t161 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t161 + Ifges(3,1) * t156 / 0.2e1) * t156) * t169 + (-t141 * t156 + t142 * t161) * mrSges(3,3)) * t169 + (t141 * mrSges(3,1) - t142 * mrSges(3,2) + Ifges(3,3) * t147 / 0.2e1 + (Ifges(3,5) * t156 + Ifges(3,6) * t161) * t169) * t147 + (-t125 * mrSges(4,1) + t124 * mrSges(4,3) + Ifges(4,4) * t136 + Ifges(4,6) * t140 + Ifges(4,2) * t135 / 0.2e1) * t135 + (-t118 * mrSges(5,1) + t110 * mrSges(5,3) + Ifges(5,4) * t127 + Ifges(5,6) * t133 + Ifges(5,2) * t126 / 0.2e1) * t126 + (-t111 * mrSges(6,1) + t103 * mrSges(6,3) + Ifges(6,4) * t122 + Ifges(6,6) * t132 + Ifges(6,2) * t121 / 0.2e1) * t121 + (-t100 * mrSges(7,1) + t99 * mrSges(7,3) + Ifges(7,4) * t113 + Ifges(7,6) * t120 + Ifges(7,2) * t112 / 0.2e1) * t112;
T  = t1;
