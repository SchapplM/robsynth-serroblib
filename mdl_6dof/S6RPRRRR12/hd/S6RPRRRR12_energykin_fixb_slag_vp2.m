% Calculate kinetic energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:20
% EndTime: 2019-03-09 07:50:21
% DurationCPUTime: 0.95s
% Computational Cost: add. (2798->129), mult. (9156->217), div. (0->0), fcn. (7926->16), ass. (0->66)
t162 = sin(pkin(14));
t166 = cos(pkin(14));
t165 = sin(pkin(6));
t182 = qJD(1) * t165;
t180 = qJ(2) * t182;
t169 = cos(pkin(6));
t181 = pkin(1) * qJD(1) * t169;
t157 = t162 * t181 + t166 * t180;
t164 = sin(pkin(7));
t168 = cos(pkin(7));
t185 = t165 * t166;
t147 = (t164 * t169 + t168 * t185) * qJD(1) * pkin(10) + t157;
t160 = t166 * t181;
t148 = t160 + (pkin(2) * t169 + (-pkin(10) * t168 - qJ(2)) * t165 * t162) * qJD(1);
t154 = qJD(2) + (-pkin(10) * t162 * t164 - pkin(2) * t166 - pkin(1)) * t182;
t173 = sin(qJ(3));
t177 = cos(qJ(3));
t183 = t168 * t177;
t186 = t164 * t177;
t137 = -t147 * t173 + t148 * t183 + t154 * t186;
t155 = qJD(3) + (-t164 * t185 + t168 * t169) * qJD(1);
t167 = cos(pkin(8));
t184 = t168 * t173;
t187 = t164 * t173;
t150 = (t169 * t187 + (t162 * t177 + t166 * t184) * t165) * qJD(1);
t191 = pkin(11) * t150;
t130 = pkin(3) * t155 - t167 * t191 + t137;
t142 = -t148 * t164 + t168 * t154;
t149 = (t169 * t186 + (-t162 * t173 + t166 * t183) * t165) * qJD(1);
t163 = sin(pkin(8));
t136 = -pkin(3) * t149 - t163 * t191 + t142;
t192 = t130 * t167 + t136 * t163;
t138 = t177 * t147 + t148 * t184 + t154 * t187;
t178 = t149 * t167 + t155 * t163;
t129 = t178 * pkin(11) + t138;
t172 = sin(qJ(4));
t176 = cos(qJ(4));
t121 = -t172 * t129 + t192 * t176;
t140 = -t150 * t172 + t178 * t176;
t122 = t176 * t129 + t192 * t172;
t143 = -t149 * t163 + t155 * t167 + qJD(4);
t118 = pkin(12) * t143 + t122;
t123 = -t130 * t163 + t167 * t136;
t141 = t150 * t176 + t178 * t172;
t120 = -pkin(4) * t140 - pkin(12) * t141 + t123;
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t114 = t175 * t118 + t171 * t120;
t113 = -t118 * t171 + t120 * t175;
t132 = -t141 * t171 + t143 * t175;
t117 = -pkin(4) * t143 - t121;
t174 = cos(qJ(6));
t170 = sin(qJ(6));
t161 = -pkin(1) * t182 + qJD(2);
t156 = -t162 * t180 + t160;
t139 = qJD(5) - t140;
t133 = t141 * t175 + t143 * t171;
t131 = qJD(6) - t132;
t125 = t133 * t174 + t139 * t170;
t124 = -t133 * t170 + t139 * t174;
t115 = -pkin(5) * t132 - pkin(13) * t133 + t117;
t112 = pkin(13) * t139 + t114;
t111 = -pkin(5) * t139 - t113;
t110 = t112 * t174 + t115 * t170;
t109 = -t112 * t170 + t115 * t174;
t1 = m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t137 ^ 2 + t138 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t156 ^ 2 + t157 ^ 2 + t161 ^ 2) / 0.2e1 + (t137 * mrSges(4,1) - t138 * mrSges(4,2) + Ifges(4,3) * t155 / 0.2e1) * t155 + (t121 * mrSges(5,1) - t122 * mrSges(5,2) + Ifges(5,3) * t143 / 0.2e1) * t143 + (t113 * mrSges(6,1) - t114 * mrSges(6,2) + Ifges(6,3) * t139 / 0.2e1) * t139 + (t109 * mrSges(7,1) - t110 * mrSges(7,2) + Ifges(7,3) * t131 / 0.2e1) * t131 + (t142 * mrSges(4,2) - t137 * mrSges(4,3) + Ifges(4,5) * t155 + Ifges(4,1) * t150 / 0.2e1) * t150 + (t123 * mrSges(5,2) - t121 * mrSges(5,3) + Ifges(5,5) * t143 + Ifges(5,1) * t141 / 0.2e1) * t141 + (t117 * mrSges(6,2) - t113 * mrSges(6,3) + Ifges(6,5) * t139 + Ifges(6,1) * t133 / 0.2e1) * t133 + (t111 * mrSges(7,2) - t109 * mrSges(7,3) + Ifges(7,5) * t131 + Ifges(7,1) * t125 / 0.2e1) * t125 + (-t142 * mrSges(4,1) + t138 * mrSges(4,3) + Ifges(4,4) * t150 + Ifges(4,6) * t155 + Ifges(4,2) * t149 / 0.2e1) * t149 + (-t123 * mrSges(5,1) + t122 * mrSges(5,3) + Ifges(5,4) * t141 + Ifges(5,6) * t143 + Ifges(5,2) * t140 / 0.2e1) * t140 + (-t117 * mrSges(6,1) + t114 * mrSges(6,3) + Ifges(6,4) * t133 + Ifges(6,6) * t139 + Ifges(6,2) * t132 / 0.2e1) * t132 + (-t111 * mrSges(7,1) + t110 * mrSges(7,3) + Ifges(7,4) * t125 + Ifges(7,6) * t131 + Ifges(7,2) * t124 / 0.2e1) * t124 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t161 * (-mrSges(3,1) * t166 + mrSges(3,2) * t162) + (Ifges(3,2) * t166 ^ 2 / 0.2e1 + (Ifges(3,4) * t166 + Ifges(3,1) * t162 / 0.2e1) * t162) * t182 + (-t156 * t162 + t157 * t166) * mrSges(3,3)) * t165 + (t156 * mrSges(3,1) - t157 * mrSges(3,2) + (Ifges(3,3) * t169 / 0.2e1 + (Ifges(3,5) * t162 + Ifges(3,6) * t166) * t165) * qJD(1)) * t169) * qJD(1);
T  = t1;
