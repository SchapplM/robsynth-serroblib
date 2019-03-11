% Calculate kinetic energy for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:22
% EndTime: 2019-03-09 04:11:22
% DurationCPUTime: 0.72s
% Computational Cost: add. (1776->124), mult. (5642->205), div. (0->0), fcn. (4776->14), ass. (0->59)
t151 = sin(pkin(12));
t155 = cos(pkin(12));
t153 = sin(pkin(6));
t169 = qJD(1) * t153;
t165 = t155 * t169;
t157 = cos(pkin(6));
t168 = qJD(1) * t157;
t167 = pkin(1) * t168;
t143 = qJ(2) * t165 + t151 * t167;
t152 = sin(pkin(7));
t156 = cos(pkin(7));
t171 = t153 * t155;
t133 = (t152 * t157 + t156 * t171) * qJD(1) * pkin(9) + t143;
t148 = t155 * t167;
t135 = t148 + (pkin(2) * t157 + (-pkin(9) * t156 - qJ(2)) * t153 * t151) * qJD(1);
t140 = qJD(2) + (-pkin(9) * t151 * t152 - pkin(2) * t155 - pkin(1)) * t169;
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t125 = -t160 * t133 + (t135 * t156 + t140 * t152) * t163;
t172 = t152 * t160;
t170 = t156 * t160;
t127 = -t135 * t152 + t156 * t140;
t166 = t151 * t169;
t136 = t160 * t166 + (-t152 * t168 - t156 * t165) * t163;
t137 = (t157 * t172 + (t151 * t163 + t155 * t170) * t153) * qJD(1);
t119 = pkin(3) * t136 - qJ(4) * t137 + t127;
t126 = t163 * t133 + t135 * t170 + t140 * t172;
t141 = qJD(3) + (-t152 * t171 + t156 * t157) * qJD(1);
t121 = qJ(4) * t141 + t126;
t150 = sin(pkin(13));
t154 = cos(pkin(13));
t111 = t154 * t119 - t121 * t150;
t129 = t137 * t154 + t141 * t150;
t108 = pkin(4) * t136 - pkin(10) * t129 + t111;
t112 = t150 * t119 + t154 * t121;
t128 = -t137 * t150 + t141 * t154;
t110 = pkin(10) * t128 + t112;
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t105 = t159 * t108 + t162 * t110;
t104 = t108 * t162 - t110 * t159;
t123 = t128 * t162 - t129 * t159;
t120 = -t141 * pkin(3) + qJD(4) - t125;
t113 = -t128 * pkin(4) + t120;
t161 = cos(qJ(6));
t158 = sin(qJ(6));
t149 = -pkin(1) * t169 + qJD(2);
t142 = -qJ(2) * t166 + t148;
t134 = qJD(5) + t136;
t124 = t128 * t159 + t129 * t162;
t122 = qJD(6) - t123;
t115 = t124 * t161 + t134 * t158;
t114 = -t124 * t158 + t134 * t161;
t106 = -t123 * pkin(5) - t124 * pkin(11) + t113;
t103 = pkin(11) * t134 + t105;
t102 = -pkin(5) * t134 - t104;
t101 = t103 * t161 + t106 * t158;
t100 = -t103 * t158 + t106 * t161;
t1 = m(3) * (t142 ^ 2 + t143 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t120 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + (t125 * mrSges(4,1) - t126 * mrSges(4,2) + Ifges(4,3) * t141 / 0.2e1) * t141 + (t104 * mrSges(6,1) - t105 * mrSges(6,2) + Ifges(6,3) * t134 / 0.2e1) * t134 + (t120 * mrSges(5,2) - t111 * mrSges(5,3) + Ifges(5,1) * t129 / 0.2e1) * t129 + (t100 * mrSges(7,1) - t101 * mrSges(7,2) + Ifges(7,3) * t122 / 0.2e1) * t122 + (t127 * mrSges(4,2) - t125 * mrSges(4,3) + Ifges(4,5) * t141 + Ifges(4,1) * t137 / 0.2e1) * t137 + (-t120 * mrSges(5,1) + t112 * mrSges(5,3) + Ifges(5,4) * t129 + Ifges(5,2) * t128 / 0.2e1) * t128 + (t113 * mrSges(6,2) - t104 * mrSges(6,3) + Ifges(6,5) * t134 + Ifges(6,1) * t124 / 0.2e1) * t124 + (t102 * mrSges(7,2) - t100 * mrSges(7,3) + Ifges(7,5) * t122 + Ifges(7,1) * t115 / 0.2e1) * t115 + (-t113 * mrSges(6,1) + t105 * mrSges(6,3) + Ifges(6,4) * t124 + Ifges(6,6) * t134 + Ifges(6,2) * t123 / 0.2e1) * t123 + (-t102 * mrSges(7,1) + t101 * mrSges(7,3) + Ifges(7,4) * t115 + Ifges(7,6) * t122 + Ifges(7,2) * t114 / 0.2e1) * t114 + (t127 * mrSges(4,1) + t111 * mrSges(5,1) - t112 * mrSges(5,2) - t126 * mrSges(4,3) - Ifges(4,4) * t137 + Ifges(5,5) * t129 - Ifges(4,6) * t141 + Ifges(5,6) * t128 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t136) * t136 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t149 * (-mrSges(3,1) * t155 + mrSges(3,2) * t151) + (Ifges(3,2) * t155 ^ 2 / 0.2e1 + (Ifges(3,4) * t155 + Ifges(3,1) * t151 / 0.2e1) * t151) * t169 + (-t142 * t151 + t143 * t155) * mrSges(3,3)) * t153 + (t142 * mrSges(3,1) - t143 * mrSges(3,2) + (Ifges(3,3) * t157 / 0.2e1 + (Ifges(3,5) * t151 + Ifges(3,6) * t155) * t153) * qJD(1)) * t157) * qJD(1);
T  = t1;
