% Calculate kinetic energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:30
% EndTime: 2019-03-08 20:50:31
% DurationCPUTime: 0.65s
% Computational Cost: add. (1076->109), mult. (3205->192), div. (0->0), fcn. (2705->16), ass. (0->57)
t167 = cos(qJ(2));
t155 = sin(pkin(6));
t176 = qJD(1) * t155;
t147 = qJD(2) * pkin(2) + t167 * t176;
t154 = sin(pkin(7));
t158 = cos(pkin(7));
t159 = cos(pkin(6));
t175 = qJD(1) * t159;
t184 = t147 * t158 + t154 * t175;
t152 = sin(pkin(14));
t162 = sin(qJ(4));
t166 = cos(qJ(4));
t156 = cos(pkin(14));
t157 = cos(pkin(8));
t178 = t156 * t157;
t153 = sin(pkin(8));
t180 = t153 * t158;
t139 = (t154 * (-t152 * t162 + t166 * t178) + t166 * t180) * qJD(2);
t163 = sin(qJ(2));
t174 = qJD(2) * t154;
t146 = qJ(3) * t174 + t163 * t176;
t136 = t156 * t146 + t184 * t152;
t128 = (t154 * t178 + t180) * qJD(2) * pkin(10) + t136;
t135 = -t146 * t152 + t184 * t156;
t182 = pkin(10) * t152;
t129 = (-t154 * t157 * t182 + pkin(3) * t158) * qJD(2) + t135;
t173 = t158 * t175 + qJD(3);
t137 = (-t147 + (-pkin(3) * t156 - t153 * t182) * qJD(2)) * t154 + t173;
t119 = -t162 * t128 + (t129 * t157 + t137 * t153) * t166;
t179 = t153 * t162;
t177 = t157 * t162;
t120 = t166 * t128 + t129 * t177 + t137 * t179;
t142 = qJD(4) + (-t153 * t154 * t156 + t157 * t158) * qJD(2);
t117 = pkin(11) * t142 + t120;
t122 = -t129 * t153 + t157 * t137;
t140 = (t158 * t179 + (t152 * t166 + t156 * t177) * t154) * qJD(2);
t121 = -pkin(4) * t139 - pkin(11) * t140 + t122;
t161 = sin(qJ(5));
t165 = cos(qJ(5));
t113 = t165 * t117 + t161 * t121;
t112 = -t117 * t161 + t121 * t165;
t131 = -t140 * t161 + t142 * t165;
t116 = -pkin(4) * t142 - t119;
t164 = cos(qJ(6));
t160 = sin(qJ(6));
t141 = -t147 * t154 + t173;
t138 = qJD(5) - t139;
t132 = t140 * t165 + t142 * t161;
t130 = qJD(6) - t131;
t124 = t132 * t164 + t138 * t160;
t123 = -t132 * t160 + t138 * t164;
t114 = -pkin(5) * t131 - pkin(12) * t132 + t116;
t111 = pkin(12) * t138 + t113;
t110 = -pkin(5) * t138 - t112;
t109 = t111 * t164 + t114 * t160;
t108 = -t111 * t160 + t114 * t164;
t1 = m(4) * (t135 ^ 2 + t136 ^ 2 + t141 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + (t119 * mrSges(5,1) - t120 * mrSges(5,2) + Ifges(5,3) * t142 / 0.2e1) * t142 + (t112 * mrSges(6,1) - t113 * mrSges(6,2) + Ifges(6,3) * t138 / 0.2e1) * t138 + (t108 * mrSges(7,1) - t109 * mrSges(7,2) + Ifges(7,3) * t130 / 0.2e1) * t130 + (t122 * mrSges(5,2) - t119 * mrSges(5,3) + Ifges(5,5) * t142 + Ifges(5,1) * t140 / 0.2e1) * t140 + (t116 * mrSges(6,2) - t112 * mrSges(6,3) + Ifges(6,5) * t138 + Ifges(6,1) * t132 / 0.2e1) * t132 + (t110 * mrSges(7,2) - t108 * mrSges(7,3) + Ifges(7,5) * t130 + Ifges(7,1) * t124 / 0.2e1) * t124 + (m(3) * (t159 ^ 2 + (t163 ^ 2 + t167 ^ 2) * t155 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t122 * mrSges(5,1) + t120 * mrSges(5,3) + Ifges(5,4) * t140 + Ifges(5,6) * t142 + Ifges(5,2) * t139 / 0.2e1) * t139 + (-t116 * mrSges(6,1) + t113 * mrSges(6,3) + Ifges(6,4) * t132 + Ifges(6,6) * t138 + Ifges(6,2) * t131 / 0.2e1) * t131 + (-t110 * mrSges(7,1) + t109 * mrSges(7,3) + Ifges(7,4) * t124 + Ifges(7,6) * t130 + Ifges(7,2) * t123 / 0.2e1) * t123 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t141 * (-mrSges(4,1) * t156 + mrSges(4,2) * t152) + (Ifges(4,2) * t156 ^ 2 / 0.2e1 + (Ifges(4,4) * t156 + Ifges(4,1) * t152 / 0.2e1) * t152) * t174 + (-t135 * t152 + t136 * t156) * mrSges(4,3)) * t154 + (t135 * mrSges(4,1) - t136 * mrSges(4,2) + (Ifges(4,3) * t158 / 0.2e1 + (Ifges(4,5) * t152 + Ifges(4,6) * t156) * t154) * qJD(2)) * t158 + (mrSges(3,1) * t167 - mrSges(3,2) * t163) * t176) * qJD(2);
T  = t1;
