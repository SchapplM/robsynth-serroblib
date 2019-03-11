% Calculate kinetic energy for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:12
% EndTime: 2019-03-09 07:38:13
% DurationCPUTime: 0.82s
% Computational Cost: add. (1790->124), mult. (5576->207), div. (0->0), fcn. (4718->14), ass. (0->59)
t151 = sin(pkin(13));
t153 = sin(pkin(6));
t160 = sin(qJ(3));
t164 = cos(qJ(3));
t154 = cos(pkin(13));
t155 = cos(pkin(7));
t172 = t154 * t155;
t152 = sin(pkin(7));
t156 = cos(pkin(6));
t174 = t152 * t156;
t139 = ((-t151 * t160 + t164 * t172) * t153 + t164 * t174) * qJD(1);
t170 = qJD(1) * t153;
t167 = qJ(2) * t170;
t169 = pkin(1) * qJD(1) * t156;
t146 = t151 * t169 + t154 * t167;
t135 = (t153 * t172 + t174) * qJD(1) * pkin(9) + t146;
t149 = t154 * t169;
t138 = t149 + (pkin(2) * t156 + (-pkin(9) * t155 - qJ(2)) * t153 * t151) * qJD(1);
t143 = qJD(2) + (-pkin(9) * t151 * t152 - pkin(2) * t154 - pkin(1)) * t170;
t123 = -t160 * t135 + (t138 * t155 + t143 * t152) * t164;
t173 = t152 * t160;
t171 = t155 * t160;
t127 = -t138 * t152 + t155 * t143;
t140 = (t156 * t173 + (t151 * t164 + t154 * t171) * t153) * qJD(1);
t118 = -pkin(3) * t139 - pkin(10) * t140 + t127;
t124 = t164 * t135 + t138 * t171 + t143 * t173;
t144 = qJD(3) + (-t152 * t153 * t154 + t155 * t156) * qJD(1);
t122 = pkin(10) * t144 + t124;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t111 = t159 * t118 + t163 * t122;
t137 = qJD(4) - t139;
t109 = pkin(11) * t137 + t111;
t121 = -t144 * pkin(3) - t123;
t130 = -t159 * t140 + t144 * t163;
t131 = t140 * t163 + t144 * t159;
t114 = -t130 * pkin(4) - t131 * pkin(11) + t121;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t105 = t162 * t109 + t158 * t114;
t104 = -t109 * t158 + t162 * t114;
t110 = t118 * t163 - t159 * t122;
t129 = qJD(5) - t130;
t108 = -pkin(4) * t137 - t110;
t161 = cos(qJ(6));
t157 = sin(qJ(6));
t150 = -pkin(1) * t170 + qJD(2);
t145 = -t151 * t167 + t149;
t128 = qJD(6) + t129;
t126 = t131 * t162 + t137 * t158;
t125 = -t131 * t158 + t137 * t162;
t116 = t125 * t157 + t126 * t161;
t115 = t125 * t161 - t126 * t157;
t106 = -pkin(5) * t125 + t108;
t103 = pkin(12) * t125 + t105;
t102 = pkin(5) * t129 - pkin(12) * t126 + t104;
t101 = t102 * t157 + t103 * t161;
t100 = t102 * t161 - t103 * t157;
t1 = m(4) * (t123 ^ 2 + t124 ^ 2 + t127 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t145 ^ 2 + t146 ^ 2 + t150 ^ 2) / 0.2e1 + (t123 * mrSges(4,1) - t124 * mrSges(4,2) + Ifges(4,3) * t144 / 0.2e1) * t144 + (t110 * mrSges(5,1) - t111 * mrSges(5,2) + Ifges(5,3) * t137 / 0.2e1) * t137 + (t104 * mrSges(6,1) - t105 * mrSges(6,2) + Ifges(6,3) * t129 / 0.2e1) * t129 + (t100 * mrSges(7,1) - t101 * mrSges(7,2) + Ifges(7,3) * t128 / 0.2e1) * t128 + (t127 * mrSges(4,2) - t123 * mrSges(4,3) + Ifges(4,5) * t144 + Ifges(4,1) * t140 / 0.2e1) * t140 + (t121 * mrSges(5,2) - t110 * mrSges(5,3) + Ifges(5,5) * t137 + Ifges(5,1) * t131 / 0.2e1) * t131 + (t108 * mrSges(6,2) - t104 * mrSges(6,3) + Ifges(6,5) * t129 + Ifges(6,1) * t126 / 0.2e1) * t126 + (t106 * mrSges(7,2) - t100 * mrSges(7,3) + Ifges(7,5) * t128 + Ifges(7,1) * t116 / 0.2e1) * t116 + (-t127 * mrSges(4,1) + t124 * mrSges(4,3) + Ifges(4,4) * t140 + Ifges(4,6) * t144 + Ifges(4,2) * t139 / 0.2e1) * t139 + (-t121 * mrSges(5,1) + t111 * mrSges(5,3) + Ifges(5,4) * t131 + Ifges(5,6) * t137 + Ifges(5,2) * t130 / 0.2e1) * t130 + (-t108 * mrSges(6,1) + t105 * mrSges(6,3) + Ifges(6,4) * t126 + Ifges(6,6) * t129 + Ifges(6,2) * t125 / 0.2e1) * t125 + (-t106 * mrSges(7,1) + t101 * mrSges(7,3) + Ifges(7,4) * t116 + Ifges(7,6) * t128 + Ifges(7,2) * t115 / 0.2e1) * t115 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t150 * (-mrSges(3,1) * t154 + mrSges(3,2) * t151) + (Ifges(3,2) * t154 ^ 2 / 0.2e1 + (Ifges(3,4) * t154 + Ifges(3,1) * t151 / 0.2e1) * t151) * t170 + (-t145 * t151 + t146 * t154) * mrSges(3,3)) * t153 + (t145 * mrSges(3,1) - t146 * mrSges(3,2) + (Ifges(3,3) * t156 / 0.2e1 + (Ifges(3,5) * t151 + Ifges(3,6) * t154) * t153) * qJD(1)) * t156) * qJD(1);
T  = t1;
