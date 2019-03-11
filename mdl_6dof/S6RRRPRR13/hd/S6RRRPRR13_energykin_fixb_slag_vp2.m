% Calculate kinetic energy for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:51
% EndTime: 2019-03-09 19:51:51
% DurationCPUTime: 0.81s
% Computational Cost: add. (2182->122), mult. (5644->197), div. (0->0), fcn. (4776->14), ass. (0->59)
t170 = cos(pkin(6)) * qJD(1);
t149 = qJD(2) + t170;
t151 = sin(pkin(7));
t154 = cos(pkin(7));
t163 = cos(qJ(2));
t152 = sin(pkin(6));
t171 = qJD(1) * t152;
t167 = t163 * t171;
t175 = t149 * t151 + t154 * t167;
t159 = sin(qJ(2));
t169 = pkin(1) * t170;
t143 = pkin(9) * t167 + t159 * t169;
t133 = t175 * pkin(10) + t143;
t148 = t163 * t169;
t168 = t159 * t171;
t135 = pkin(2) * t149 + t148 + (-pkin(10) * t154 - pkin(9)) * t168;
t140 = (-pkin(10) * t151 * t159 - pkin(2) * t163 - pkin(1)) * t171;
t158 = sin(qJ(3));
t162 = cos(qJ(3));
t125 = -t158 * t133 + (t135 * t154 + t140 * t151) * t162;
t173 = t151 * t158;
t172 = t154 * t158;
t127 = -t135 * t151 + t154 * t140;
t136 = t158 * t168 - t175 * t162;
t137 = t149 * t173 + (t159 * t162 + t163 * t172) * t171;
t118 = pkin(3) * t136 - qJ(4) * t137 + t127;
t126 = t162 * t133 + t135 * t172 + t140 * t173;
t141 = t149 * t154 - t151 * t167 + qJD(3);
t121 = qJ(4) * t141 + t126;
t150 = sin(pkin(13));
t153 = cos(pkin(13));
t111 = t153 * t118 - t121 * t150;
t129 = t137 * t153 + t141 * t150;
t108 = pkin(4) * t136 - pkin(11) * t129 + t111;
t112 = t150 * t118 + t153 * t121;
t128 = -t137 * t150 + t141 * t153;
t110 = pkin(11) * t128 + t112;
t157 = sin(qJ(5));
t161 = cos(qJ(5));
t105 = t157 * t108 + t161 * t110;
t104 = t108 * t161 - t110 * t157;
t123 = t128 * t161 - t129 * t157;
t120 = -pkin(3) * t141 + qJD(4) - t125;
t113 = -pkin(4) * t128 + t120;
t164 = qJD(1) ^ 2;
t160 = cos(qJ(6));
t156 = sin(qJ(6));
t142 = -pkin(9) * t168 + t148;
t134 = qJD(5) + t136;
t124 = t128 * t157 + t129 * t161;
t122 = qJD(6) - t123;
t115 = t124 * t160 + t134 * t156;
t114 = -t124 * t156 + t134 * t160;
t106 = -pkin(5) * t123 - pkin(12) * t124 + t113;
t103 = pkin(12) * t134 + t105;
t102 = -pkin(5) * t134 - t104;
t101 = t103 * t160 + t106 * t156;
t100 = -t103 * t156 + t106 * t160;
t1 = m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + t164 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t152 ^ 2 * t164 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + (t125 * mrSges(4,1) - t126 * mrSges(4,2) + Ifges(4,3) * t141 / 0.2e1) * t141 + (t104 * mrSges(6,1) - t105 * mrSges(6,2) + Ifges(6,3) * t134 / 0.2e1) * t134 + (t120 * mrSges(5,2) - t111 * mrSges(5,3) + Ifges(5,1) * t129 / 0.2e1) * t129 + (t100 * mrSges(7,1) - t101 * mrSges(7,2) + Ifges(7,3) * t122 / 0.2e1) * t122 + (t127 * mrSges(4,2) - t125 * mrSges(4,3) + Ifges(4,5) * t141 + Ifges(4,1) * t137 / 0.2e1) * t137 + (-t120 * mrSges(5,1) + t112 * mrSges(5,3) + Ifges(5,4) * t129 + Ifges(5,2) * t128 / 0.2e1) * t128 + (t113 * mrSges(6,2) - t104 * mrSges(6,3) + Ifges(6,5) * t134 + Ifges(6,1) * t124 / 0.2e1) * t124 + (t102 * mrSges(7,2) - t100 * mrSges(7,3) + Ifges(7,5) * t122 + Ifges(7,1) * t115 / 0.2e1) * t115 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t163 / 0.2e1) * t163 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t163 + Ifges(3,1) * t159 / 0.2e1) * t159) * t171 + (-t142 * t159 + t143 * t163) * mrSges(3,3)) * t171 + (t142 * mrSges(3,1) - t143 * mrSges(3,2) + Ifges(3,3) * t149 / 0.2e1 + (Ifges(3,5) * t159 + Ifges(3,6) * t163) * t171) * t149 + (-t113 * mrSges(6,1) + t105 * mrSges(6,3) + Ifges(6,4) * t124 + Ifges(6,6) * t134 + Ifges(6,2) * t123 / 0.2e1) * t123 + (-t102 * mrSges(7,1) + t101 * mrSges(7,3) + Ifges(7,4) * t115 + Ifges(7,6) * t122 + Ifges(7,2) * t114 / 0.2e1) * t114 + (t127 * mrSges(4,1) + t111 * mrSges(5,1) - t112 * mrSges(5,2) - t126 * mrSges(4,3) - Ifges(4,4) * t137 + Ifges(5,5) * t129 - Ifges(4,6) * t141 + Ifges(5,6) * t128 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t136) * t136;
T  = t1;
