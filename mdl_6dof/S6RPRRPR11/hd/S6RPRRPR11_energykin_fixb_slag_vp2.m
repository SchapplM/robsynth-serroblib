% Calculate kinetic energy for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:00
% EndTime: 2019-03-09 05:39:01
% DurationCPUTime: 0.78s
% Computational Cost: add. (1778->124), mult. (5576->205), div. (0->0), fcn. (4718->14), ass. (0->58)
t150 = sin(pkin(12));
t152 = sin(pkin(6));
t159 = sin(qJ(3));
t161 = cos(qJ(3));
t154 = cos(pkin(12));
t155 = cos(pkin(7));
t169 = t154 * t155;
t151 = sin(pkin(7));
t156 = cos(pkin(6));
t171 = t151 * t156;
t136 = (t152 * (-t150 * t159 + t161 * t169) + t161 * t171) * qJD(1);
t167 = qJD(1) * t152;
t164 = qJ(2) * t167;
t166 = pkin(1) * qJD(1) * t156;
t144 = t150 * t166 + t154 * t164;
t133 = (t152 * t169 + t171) * qJD(1) * pkin(9) + t144;
t147 = t154 * t166;
t135 = t147 + (pkin(2) * t156 + (-pkin(9) * t155 - qJ(2)) * t152 * t150) * qJD(1);
t140 = qJD(2) + (-pkin(9) * t150 * t151 - pkin(2) * t154 - pkin(1)) * t167;
t122 = -t159 * t133 + (t135 * t155 + t140 * t151) * t161;
t172 = cos(qJ(4));
t170 = t151 * t159;
t168 = t155 * t159;
t126 = -t135 * t151 + t155 * t140;
t137 = (t156 * t170 + (t150 * t161 + t154 * t168) * t152) * qJD(1);
t117 = -pkin(3) * t136 - pkin(10) * t137 + t126;
t123 = t161 * t133 + t135 * t168 + t140 * t170;
t142 = qJD(3) + (-t151 * t152 * t154 + t155 * t156) * qJD(1);
t121 = pkin(10) * t142 + t123;
t158 = sin(qJ(4));
t110 = t158 * t117 + t172 * t121;
t134 = qJD(4) - t136;
t108 = qJ(5) * t134 + t110;
t120 = -t142 * pkin(3) - t122;
t128 = t137 * t158 - t172 * t142;
t129 = t137 * t172 + t158 * t142;
t113 = t128 * pkin(4) - t129 * qJ(5) + t120;
t149 = sin(pkin(13));
t153 = cos(pkin(13));
t104 = t153 * t108 + t149 * t113;
t103 = -t108 * t149 + t153 * t113;
t109 = t117 * t172 - t158 * t121;
t107 = -t134 * pkin(4) + qJD(5) - t109;
t160 = cos(qJ(6));
t157 = sin(qJ(6));
t148 = -pkin(1) * t167 + qJD(2);
t143 = -t150 * t164 + t147;
t127 = qJD(6) + t128;
t125 = t129 * t153 + t134 * t149;
t124 = -t129 * t149 + t134 * t153;
t115 = t124 * t157 + t125 * t160;
t114 = t124 * t160 - t125 * t157;
t105 = -t124 * pkin(5) + t107;
t102 = pkin(11) * t124 + t104;
t101 = pkin(5) * t128 - pkin(11) * t125 + t103;
t100 = t101 * t157 + t102 * t160;
t99 = t101 * t160 - t102 * t157;
t1 = m(5) * (t109 ^ 2 + t110 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t143 ^ 2 + t144 ^ 2 + t148 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t126 ^ 2) / 0.2e1 + (t122 * mrSges(4,1) - t123 * mrSges(4,2) + Ifges(4,3) * t142 / 0.2e1) * t142 + (t109 * mrSges(5,1) - t110 * mrSges(5,2) + Ifges(5,3) * t134 / 0.2e1) * t134 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t127 / 0.2e1) * t127 + (t107 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,1) * t125 / 0.2e1) * t125 + (t126 * mrSges(4,2) - t122 * mrSges(4,3) + Ifges(4,5) * t142 + Ifges(4,1) * t137 / 0.2e1) * t137 + (t120 * mrSges(5,2) - t109 * mrSges(5,3) + Ifges(5,5) * t134 + Ifges(5,1) * t129 / 0.2e1) * t129 + (-t107 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t125 + Ifges(6,2) * t124 / 0.2e1) * t124 + (t105 * mrSges(7,2) - t99 * mrSges(7,3) + Ifges(7,5) * t127 + Ifges(7,1) * t115 / 0.2e1) * t115 + (-t126 * mrSges(4,1) + t123 * mrSges(4,3) + Ifges(4,4) * t137 + Ifges(4,6) * t142 + Ifges(4,2) * t136 / 0.2e1) * t136 + (-t105 * mrSges(7,1) + t100 * mrSges(7,3) + Ifges(7,4) * t115 + Ifges(7,6) * t127 + Ifges(7,2) * t114 / 0.2e1) * t114 + (t120 * mrSges(5,1) + t103 * mrSges(6,1) - t104 * mrSges(6,2) - t110 * mrSges(5,3) - Ifges(5,4) * t129 + Ifges(6,5) * t125 - Ifges(5,6) * t134 + Ifges(6,6) * t124 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t128) * t128 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t148 * (-mrSges(3,1) * t154 + mrSges(3,2) * t150) + (Ifges(3,2) * t154 ^ 2 / 0.2e1 + (Ifges(3,4) * t154 + Ifges(3,1) * t150 / 0.2e1) * t150) * t167 + (-t143 * t150 + t144 * t154) * mrSges(3,3)) * t152 + (t143 * mrSges(3,1) - t144 * mrSges(3,2) + (Ifges(3,3) * t156 / 0.2e1 + (Ifges(3,5) * t150 + Ifges(3,6) * t154) * t152) * qJD(1)) * t156) * qJD(1);
T  = t1;
