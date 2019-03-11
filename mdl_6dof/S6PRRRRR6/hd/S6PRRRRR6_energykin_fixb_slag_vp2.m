% Calculate kinetic energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:14:00
% EndTime: 2019-03-09 01:14:00
% DurationCPUTime: 0.65s
% Computational Cost: add. (1254->110), mult. (3205->187), div. (0->0), fcn. (2705->16), ass. (0->59)
t168 = cos(qJ(2));
t155 = sin(pkin(6));
t177 = qJD(1) * t155;
t147 = qJD(2) * pkin(2) + t168 * t177;
t154 = sin(pkin(7));
t157 = cos(pkin(7));
t158 = cos(pkin(6));
t176 = qJD(1) * t158;
t184 = t147 * t157 + t154 * t176;
t163 = sin(qJ(2));
t175 = qJD(2) * t154;
t146 = pkin(10) * t175 + t163 * t177;
t162 = sin(qJ(3));
t167 = cos(qJ(3));
t136 = t167 * t146 + t184 * t162;
t156 = cos(pkin(8));
t172 = t167 * t175;
t151 = qJD(2) * t157 + qJD(3);
t153 = sin(pkin(8));
t181 = t151 * t153;
t129 = (t156 * t172 + t181) * pkin(11) + t136;
t150 = t157 * t176;
t137 = t150 + (-t147 + (-pkin(11) * t153 * t162 - pkin(3) * t167) * qJD(2)) * t154;
t161 = sin(qJ(4));
t166 = cos(qJ(4));
t178 = t184 * t167;
t132 = pkin(3) * t151 + (-pkin(11) * t156 * t175 - t146) * t162 + t178;
t183 = t132 * t156;
t120 = -t161 * t129 + (t137 * t153 + t183) * t166;
t179 = t156 * t167;
t139 = (-t161 * t162 + t166 * t179) * t175 + t166 * t181;
t180 = t153 * t161;
t121 = t166 * t129 + t137 * t180 + t161 * t183;
t142 = t151 * t156 - t153 * t172 + qJD(4);
t117 = pkin(12) * t142 + t121;
t122 = -t132 * t153 + t156 * t137;
t140 = t151 * t180 + (t161 * t179 + t162 * t166) * t175;
t119 = -pkin(4) * t139 - pkin(12) * t140 + t122;
t160 = sin(qJ(5));
t165 = cos(qJ(5));
t113 = t165 * t117 + t160 * t119;
t112 = -t117 * t160 + t119 * t165;
t130 = -t140 * t160 + t142 * t165;
t116 = -pkin(4) * t142 - t120;
t164 = cos(qJ(6));
t159 = sin(qJ(6));
t141 = -t147 * t154 + t150;
t138 = qJD(5) - t139;
t135 = -t146 * t162 + t178;
t131 = t140 * t165 + t142 * t160;
t128 = qJD(6) - t130;
t124 = t131 * t164 + t138 * t159;
t123 = -t131 * t159 + t138 * t164;
t114 = -pkin(5) * t130 - pkin(13) * t131 + t116;
t111 = pkin(13) * t138 + t113;
t110 = -pkin(5) * t138 - t112;
t109 = t111 * t164 + t114 * t159;
t108 = -t111 * t159 + t114 * t164;
t1 = m(4) * (t135 ^ 2 + t136 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + (t135 * mrSges(4,1) - t136 * mrSges(4,2) + Ifges(4,3) * t151 / 0.2e1) * t151 + (t120 * mrSges(5,1) - t121 * mrSges(5,2) + Ifges(5,3) * t142 / 0.2e1) * t142 + (t112 * mrSges(6,1) - t113 * mrSges(6,2) + Ifges(6,3) * t138 / 0.2e1) * t138 + (t108 * mrSges(7,1) - t109 * mrSges(7,2) + Ifges(7,3) * t128 / 0.2e1) * t128 + (t122 * mrSges(5,2) - t120 * mrSges(5,3) + Ifges(5,5) * t142 + Ifges(5,1) * t140 / 0.2e1) * t140 + (t116 * mrSges(6,2) - t112 * mrSges(6,3) + Ifges(6,5) * t138 + Ifges(6,1) * t131 / 0.2e1) * t131 + (t110 * mrSges(7,2) - t108 * mrSges(7,3) + Ifges(7,5) * t128 + Ifges(7,1) * t124 / 0.2e1) * t124 + (m(3) * (t158 ^ 2 + (t163 ^ 2 + t168 ^ 2) * t155 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t122 * mrSges(5,1) + t121 * mrSges(5,3) + Ifges(5,4) * t140 + Ifges(5,6) * t142 + Ifges(5,2) * t139 / 0.2e1) * t139 + (-t116 * mrSges(6,1) + t113 * mrSges(6,3) + Ifges(6,4) * t131 + Ifges(6,6) * t138 + Ifges(6,2) * t130 / 0.2e1) * t130 + (-t110 * mrSges(7,1) + t109 * mrSges(7,3) + Ifges(7,4) * t124 + Ifges(7,6) * t128 + Ifges(7,2) * t123 / 0.2e1) * t123 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t168 - mrSges(3,2) * t163) * t177 + (t141 * (-mrSges(4,1) * t167 + mrSges(4,2) * t162) + (Ifges(4,2) * t167 ^ 2 / 0.2e1 + (Ifges(4,4) * t167 + Ifges(4,1) * t162 / 0.2e1) * t162) * t175 + (-t135 * t162 + t136 * t167) * mrSges(4,3) + t151 * (Ifges(4,5) * t162 + Ifges(4,6) * t167)) * t154) * qJD(2);
T  = t1;
