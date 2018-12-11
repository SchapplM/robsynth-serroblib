% Calculate kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:47
% EndTime: 2018-12-10 18:09:48
% DurationCPUTime: 1.11s
% Computational Cost: add. (3268->128), mult. (9158->212), div. (0->0), fcn. (7926->16), ass. (0->67)
t173 = sin(qJ(2));
t177 = cos(qJ(2));
t165 = sin(pkin(6));
t184 = t165 * qJD(1);
t181 = t177 * t184;
t185 = cos(pkin(6)) * qJD(1);
t183 = pkin(1) * t185;
t157 = pkin(10) * t181 + t173 * t183;
t161 = qJD(2) + t185;
t164 = sin(pkin(7));
t168 = cos(pkin(7));
t147 = (t161 * t164 + t168 * t181) * qJ(3) + t157;
t154 = (-qJ(3) * t164 * t173 - pkin(2) * t177 - pkin(1)) * t184;
t162 = sin(pkin(14));
t166 = cos(pkin(14));
t187 = t164 * t166;
t160 = t177 * t183;
t182 = t173 * t184;
t150 = pkin(2) * t161 + t160 + (-qJ(3) * t168 - pkin(10)) * t182;
t189 = t150 * t168;
t137 = -t147 * t162 + t154 * t187 + t166 * t189;
t155 = t161 * t168 - t164 * t181;
t167 = cos(pkin(8));
t186 = t168 * t177;
t188 = t162 * t164;
t149 = t161 * t188 + (t162 * t186 + t166 * t173) * t184;
t193 = pkin(11) * t149;
t133 = pkin(3) * t155 - t167 * t193 + t137;
t142 = -t150 * t164 + t168 * t154 + qJD(3);
t148 = t161 * t187 + (-t162 * t173 + t166 * t186) * t184;
t163 = sin(pkin(8));
t136 = -pkin(3) * t148 - t163 * t193 + t142;
t194 = t133 * t167 + t136 * t163;
t138 = t166 * t147 + t154 * t188 + t162 * t189;
t179 = t148 * t167 + t155 * t163;
t130 = t179 * pkin(11) + t138;
t172 = sin(qJ(4));
t176 = cos(qJ(4));
t121 = -t172 * t130 + t176 * t194;
t140 = -t149 * t172 + t179 * t176;
t122 = t176 * t130 + t172 * t194;
t143 = -t148 * t163 + t155 * t167 + qJD(4);
t118 = pkin(12) * t143 + t122;
t123 = -t133 * t163 + t167 * t136;
t141 = t149 * t176 + t179 * t172;
t120 = -pkin(4) * t140 - pkin(12) * t141 + t123;
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t114 = t175 * t118 + t171 * t120;
t113 = -t118 * t171 + t120 * t175;
t131 = -t141 * t171 + t143 * t175;
t117 = -pkin(4) * t143 - t121;
t178 = qJD(1) ^ 2;
t174 = cos(qJ(6));
t170 = sin(qJ(6));
t156 = -pkin(10) * t182 + t160;
t139 = qJD(5) - t140;
t132 = t141 * t175 + t143 * t171;
t129 = qJD(6) - t131;
t125 = t132 * t174 + t139 * t170;
t124 = -t132 * t170 + t139 * t174;
t115 = -pkin(5) * t131 - pkin(13) * t132 + t117;
t112 = pkin(13) * t139 + t114;
t111 = -pkin(5) * t139 - t113;
t110 = t112 * t174 + t115 * t170;
t109 = -t112 * t170 + t115 * t174;
t1 = m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t137 ^ 2 + t138 ^ 2 + t142 ^ 2) / 0.2e1 + t178 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t165 ^ 2 * t178 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + (t137 * mrSges(4,1) - t138 * mrSges(4,2) + Ifges(4,3) * t155 / 0.2e1) * t155 + (t121 * mrSges(5,1) - t122 * mrSges(5,2) + Ifges(5,3) * t143 / 0.2e1) * t143 + (t113 * mrSges(6,1) - t114 * mrSges(6,2) + Ifges(6,3) * t139 / 0.2e1) * t139 + (t109 * mrSges(7,1) - t110 * mrSges(7,2) + Ifges(7,3) * t129 / 0.2e1) * t129 + (t142 * mrSges(4,2) - t137 * mrSges(4,3) + Ifges(4,5) * t155 + Ifges(4,1) * t149 / 0.2e1) * t149 + (t123 * mrSges(5,2) - t121 * mrSges(5,3) + Ifges(5,5) * t143 + Ifges(5,1) * t141 / 0.2e1) * t141 + (t117 * mrSges(6,2) - t113 * mrSges(6,3) + Ifges(6,5) * t139 + Ifges(6,1) * t132 / 0.2e1) * t132 + (t111 * mrSges(7,2) - t109 * mrSges(7,3) + Ifges(7,5) * t129 + Ifges(7,1) * t125 / 0.2e1) * t125 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t177 / 0.2e1) * t177 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t177 + Ifges(3,1) * t173 / 0.2e1) * t173) * t184 + (-t156 * t173 + t157 * t177) * mrSges(3,3)) * t184 + (t156 * mrSges(3,1) - t157 * mrSges(3,2) + Ifges(3,3) * t161 / 0.2e1 + (Ifges(3,5) * t173 + Ifges(3,6) * t177) * t184) * t161 + (-t142 * mrSges(4,1) + t138 * mrSges(4,3) + Ifges(4,4) * t149 + Ifges(4,6) * t155 + Ifges(4,2) * t148 / 0.2e1) * t148 + (-t123 * mrSges(5,1) + t122 * mrSges(5,3) + Ifges(5,4) * t141 + Ifges(5,6) * t143 + Ifges(5,2) * t140 / 0.2e1) * t140 + (-t117 * mrSges(6,1) + t114 * mrSges(6,3) + Ifges(6,4) * t132 + Ifges(6,6) * t139 + Ifges(6,2) * t131 / 0.2e1) * t131 + (-t111 * mrSges(7,1) + t110 * mrSges(7,3) + Ifges(7,4) * t125 + Ifges(7,6) * t129 + Ifges(7,2) * t124 / 0.2e1) * t124;
T  = t1;
