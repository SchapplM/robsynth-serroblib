% Calculate kinetic energy for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:36
% EndTime: 2019-03-10 00:08:37
% DurationCPUTime: 0.78s
% Computational Cost: add. (2176->123), mult. (5578->199), div. (0->0), fcn. (4718->14), ass. (0->60)
t158 = sin(qJ(2));
t161 = cos(qJ(2));
t151 = sin(pkin(6));
t170 = qJD(1) * t151;
t165 = t161 * t170;
t169 = cos(pkin(6)) * qJD(1);
t168 = pkin(1) * t169;
t144 = pkin(9) * t165 + t158 * t168;
t153 = cos(pkin(7));
t148 = qJD(2) + t169;
t150 = sin(pkin(7));
t173 = t148 * t150;
t133 = (t153 * t165 + t173) * pkin(10) + t144;
t141 = (-pkin(10) * t150 * t158 - pkin(2) * t161 - pkin(1)) * t170;
t157 = sin(qJ(3));
t160 = cos(qJ(3));
t147 = t161 * t168;
t166 = t158 * t170;
t135 = pkin(2) * t148 + t147 + (-pkin(10) * t153 - pkin(9)) * t166;
t174 = t135 * t153;
t122 = -t157 * t133 + (t141 * t150 + t174) * t160;
t171 = t153 * t161;
t136 = (-t157 * t158 + t160 * t171) * t170 + t160 * t173;
t175 = cos(qJ(4));
t172 = t150 * t157;
t126 = -t135 * t150 + t153 * t141;
t137 = t148 * t172 + (t157 * t171 + t158 * t160) * t170;
t117 = -pkin(3) * t136 - pkin(11) * t137 + t126;
t123 = t160 * t133 + t141 * t172 + t157 * t174;
t142 = t148 * t153 - t150 * t165 + qJD(3);
t121 = pkin(11) * t142 + t123;
t156 = sin(qJ(4));
t110 = t156 * t117 + t175 * t121;
t134 = qJD(4) - t136;
t108 = qJ(5) * t134 + t110;
t120 = -pkin(3) * t142 - t122;
t128 = t137 * t156 - t175 * t142;
t129 = t175 * t137 + t156 * t142;
t113 = pkin(4) * t128 - qJ(5) * t129 + t120;
t149 = sin(pkin(13));
t152 = cos(pkin(13));
t104 = t152 * t108 + t149 * t113;
t103 = -t108 * t149 + t152 * t113;
t109 = t175 * t117 - t156 * t121;
t107 = -t134 * pkin(4) + qJD(5) - t109;
t162 = qJD(1) ^ 2;
t159 = cos(qJ(6));
t155 = sin(qJ(6));
t143 = -pkin(9) * t166 + t147;
t127 = qJD(6) + t128;
t125 = t129 * t152 + t134 * t149;
t124 = -t129 * t149 + t134 * t152;
t115 = t124 * t155 + t125 * t159;
t114 = t124 * t159 - t125 * t155;
t105 = -t124 * pkin(5) + t107;
t102 = pkin(12) * t124 + t104;
t101 = pkin(5) * t128 - pkin(12) * t125 + t103;
t100 = t101 * t155 + t102 * t159;
t99 = t101 * t159 - t102 * t155;
t1 = m(6) * (t103 ^ 2 + t104 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + t162 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t151 ^ 2 * t162 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t126 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t120 ^ 2) / 0.2e1 + (t122 * mrSges(4,1) - t123 * mrSges(4,2) + Ifges(4,3) * t142 / 0.2e1) * t142 + (t109 * mrSges(5,1) - t110 * mrSges(5,2) + Ifges(5,3) * t134 / 0.2e1) * t134 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t127 / 0.2e1) * t127 + (t107 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,1) * t125 / 0.2e1) * t125 + (t126 * mrSges(4,2) - t122 * mrSges(4,3) + Ifges(4,5) * t142 + Ifges(4,1) * t137 / 0.2e1) * t137 + (t120 * mrSges(5,2) - t109 * mrSges(5,3) + Ifges(5,5) * t134 + Ifges(5,1) * t129 / 0.2e1) * t129 + (-t107 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t125 + Ifges(6,2) * t124 / 0.2e1) * t124 + (t105 * mrSges(7,2) - t99 * mrSges(7,3) + Ifges(7,5) * t127 + Ifges(7,1) * t115 / 0.2e1) * t115 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t161 / 0.2e1) * t161 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t161 + Ifges(3,1) * t158 / 0.2e1) * t158) * t170 + (-t143 * t158 + t144 * t161) * mrSges(3,3)) * t170 + (t143 * mrSges(3,1) - t144 * mrSges(3,2) + Ifges(3,3) * t148 / 0.2e1 + (Ifges(3,5) * t158 + Ifges(3,6) * t161) * t170) * t148 + (-t126 * mrSges(4,1) + t123 * mrSges(4,3) + Ifges(4,4) * t137 + Ifges(4,6) * t142 + Ifges(4,2) * t136 / 0.2e1) * t136 + (-t105 * mrSges(7,1) + t100 * mrSges(7,3) + Ifges(7,4) * t115 + Ifges(7,6) * t127 + Ifges(7,2) * t114 / 0.2e1) * t114 + (t120 * mrSges(5,1) + t103 * mrSges(6,1) - t104 * mrSges(6,2) - t110 * mrSges(5,3) - Ifges(5,4) * t129 + Ifges(6,5) * t125 - Ifges(5,6) * t134 + Ifges(6,6) * t124 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t128) * t128;
T  = t1;
