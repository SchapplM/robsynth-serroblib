% Calculate kinetic energy for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:38
% EndTime: 2019-03-10 05:16:39
% DurationCPUTime: 0.86s
% Computational Cost: add. (2188->123), mult. (5578->201), div. (0->0), fcn. (4718->14), ass. (0->61)
t159 = sin(qJ(2));
t164 = cos(qJ(2));
t152 = sin(pkin(6));
t172 = t152 * qJD(1);
t168 = t164 * t172;
t173 = cos(pkin(6)) * qJD(1);
t171 = pkin(1) * t173;
t146 = pkin(9) * t168 + t159 * t171;
t153 = cos(pkin(7));
t150 = qJD(2) + t173;
t151 = sin(pkin(7));
t176 = t150 * t151;
t135 = (t153 * t168 + t176) * pkin(10) + t146;
t143 = (-pkin(10) * t151 * t159 - pkin(2) * t164 - pkin(1)) * t172;
t158 = sin(qJ(3));
t163 = cos(qJ(3));
t149 = t164 * t171;
t169 = t159 * t172;
t138 = pkin(2) * t150 + t149 + (-pkin(10) * t153 - pkin(9)) * t169;
t177 = t138 * t153;
t123 = -t158 * t135 + t163 * (t143 * t151 + t177);
t174 = t153 * t164;
t139 = (-t158 * t159 + t163 * t174) * t172 + t163 * t176;
t175 = t151 * t158;
t127 = -t138 * t151 + t153 * t143;
t140 = t150 * t175 + (t158 * t174 + t159 * t163) * t172;
t118 = -pkin(3) * t139 - pkin(11) * t140 + t127;
t124 = t163 * t135 + t143 * t175 + t158 * t177;
t144 = t150 * t153 - t151 * t168 + qJD(3);
t122 = pkin(11) * t144 + t124;
t157 = sin(qJ(4));
t162 = cos(qJ(4));
t111 = t157 * t118 + t162 * t122;
t137 = qJD(4) - t139;
t109 = pkin(12) * t137 + t111;
t121 = -pkin(3) * t144 - t123;
t130 = -t157 * t140 + t144 * t162;
t131 = t140 * t162 + t144 * t157;
t114 = -pkin(4) * t130 - pkin(12) * t131 + t121;
t156 = sin(qJ(5));
t161 = cos(qJ(5));
t105 = t161 * t109 + t156 * t114;
t104 = -t109 * t156 + t161 * t114;
t110 = t118 * t162 - t157 * t122;
t129 = qJD(5) - t130;
t108 = -pkin(4) * t137 - t110;
t165 = qJD(1) ^ 2;
t160 = cos(qJ(6));
t155 = sin(qJ(6));
t145 = -pkin(9) * t169 + t149;
t128 = qJD(6) + t129;
t126 = t131 * t161 + t137 * t156;
t125 = -t131 * t156 + t137 * t161;
t116 = t125 * t155 + t126 * t160;
t115 = t125 * t160 - t126 * t155;
t106 = -pkin(5) * t125 + t108;
t103 = pkin(13) * t125 + t105;
t102 = pkin(5) * t129 - pkin(13) * t126 + t104;
t101 = t102 * t155 + t103 * t160;
t100 = t102 * t160 - t103 * t155;
t1 = m(7) * (t100 ^ 2 + t101 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t108 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t152 ^ 2 * t165 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + t165 * Ifges(2,3) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t121 ^ 2) / 0.2e1 + (t123 * mrSges(4,1) - t124 * mrSges(4,2) + Ifges(4,3) * t144 / 0.2e1) * t144 + (t110 * mrSges(5,1) - t111 * mrSges(5,2) + Ifges(5,3) * t137 / 0.2e1) * t137 + (t104 * mrSges(6,1) - t105 * mrSges(6,2) + Ifges(6,3) * t129 / 0.2e1) * t129 + (t100 * mrSges(7,1) - t101 * mrSges(7,2) + Ifges(7,3) * t128 / 0.2e1) * t128 + (t127 * mrSges(4,2) - t123 * mrSges(4,3) + Ifges(4,5) * t144 + Ifges(4,1) * t140 / 0.2e1) * t140 + (t121 * mrSges(5,2) - t110 * mrSges(5,3) + Ifges(5,5) * t137 + Ifges(5,1) * t131 / 0.2e1) * t131 + (t108 * mrSges(6,2) - t104 * mrSges(6,3) + Ifges(6,5) * t129 + Ifges(6,1) * t126 / 0.2e1) * t126 + (t106 * mrSges(7,2) - t100 * mrSges(7,3) + Ifges(7,5) * t128 + Ifges(7,1) * t116 / 0.2e1) * t116 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t164 / 0.2e1) * t164 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t164 + Ifges(3,1) * t159 / 0.2e1) * t159) * t172 + (-t145 * t159 + t146 * t164) * mrSges(3,3)) * t172 + (t145 * mrSges(3,1) - t146 * mrSges(3,2) + Ifges(3,3) * t150 / 0.2e1 + (Ifges(3,5) * t159 + Ifges(3,6) * t164) * t172) * t150 + (-t127 * mrSges(4,1) + t124 * mrSges(4,3) + Ifges(4,4) * t140 + Ifges(4,6) * t144 + Ifges(4,2) * t139 / 0.2e1) * t139 + (-t121 * mrSges(5,1) + t111 * mrSges(5,3) + Ifges(5,4) * t131 + Ifges(5,6) * t137 + Ifges(5,2) * t130 / 0.2e1) * t130 + (-t108 * mrSges(6,1) + t105 * mrSges(6,3) + Ifges(6,4) * t126 + Ifges(6,6) * t129 + Ifges(6,2) * t125 / 0.2e1) * t125 + (-t106 * mrSges(7,1) + t101 * mrSges(7,3) + Ifges(7,4) * t116 + Ifges(7,6) * t128 + Ifges(7,2) * t115 / 0.2e1) * t115;
T  = t1;
