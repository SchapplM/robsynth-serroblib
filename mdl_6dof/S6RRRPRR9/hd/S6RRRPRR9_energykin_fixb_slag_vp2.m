% Calculate kinetic energy for
% S6RRRPRR9
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:27
% EndTime: 2019-03-09 18:56:28
% DurationCPUTime: 0.73s
% Computational Cost: add. (2082->123), mult. (5674->200), div. (0->0), fcn. (4810->14), ass. (0->60)
t158 = cos(qJ(2));
t163 = cos(pkin(6)) * qJD(1);
t162 = pkin(1) * t163;
t143 = t158 * t162;
t144 = qJD(2) + t163;
t149 = cos(pkin(7));
t154 = sin(qJ(2));
t147 = sin(pkin(6));
t164 = qJD(1) * t147;
t161 = t154 * t164;
t131 = pkin(2) * t144 + t143 + (-pkin(10) * t149 - pkin(9)) * t161;
t168 = t131 * t149;
t146 = sin(pkin(7));
t153 = sin(qJ(3));
t167 = t146 * t153;
t157 = cos(qJ(3));
t166 = t146 * t157;
t165 = t149 * t158;
t160 = t158 * t164;
t140 = pkin(9) * t160 + t154 * t162;
t130 = (t144 * t146 + t149 * t160) * pkin(10) + t140;
t137 = (-pkin(10) * t146 * t154 - pkin(2) * t158 - pkin(1)) * t164;
t117 = -t130 * t153 + t137 * t166 + t157 * t168;
t133 = t144 * t167 + (t153 * t165 + t154 * t157) * t164;
t138 = t144 * t149 - t146 * t160 + qJD(3);
t113 = pkin(3) * t138 - qJ(4) * t133 + t117;
t118 = t157 * t130 + t137 * t167 + t153 * t168;
t132 = t144 * t166 + (-t153 * t154 + t157 * t165) * t164;
t116 = qJ(4) * t132 + t118;
t145 = sin(pkin(13));
t148 = cos(pkin(13));
t107 = t145 * t113 + t148 * t116;
t105 = pkin(11) * t138 + t107;
t126 = -t131 * t146 + t149 * t137;
t122 = -pkin(3) * t132 + qJD(4) + t126;
t124 = t132 * t148 - t133 * t145;
t125 = t132 * t145 + t133 * t148;
t109 = -pkin(4) * t124 - pkin(11) * t125 + t122;
t152 = sin(qJ(5));
t156 = cos(qJ(5));
t101 = t156 * t105 + t152 * t109;
t106 = t113 * t148 - t145 * t116;
t100 = -t105 * t152 + t109 * t156;
t120 = -t125 * t152 + t138 * t156;
t104 = -pkin(4) * t138 - t106;
t159 = qJD(1) ^ 2;
t155 = cos(qJ(6));
t151 = sin(qJ(6));
t139 = -pkin(9) * t161 + t143;
t123 = qJD(5) - t124;
t121 = t125 * t156 + t138 * t152;
t119 = qJD(6) - t120;
t111 = t121 * t155 + t123 * t151;
t110 = -t121 * t151 + t123 * t155;
t102 = -pkin(5) * t120 - pkin(12) * t121 + t104;
t99 = pkin(12) * t123 + t101;
t98 = -pkin(5) * t123 - t100;
t97 = t102 * t151 + t155 * t99;
t96 = t102 * t155 - t151 * t99;
t1 = m(3) * (pkin(1) ^ 2 * t147 ^ 2 * t159 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + t159 * Ifges(2,3) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t126 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t104 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + (t126 * mrSges(4,2) - t117 * mrSges(4,3) + Ifges(4,1) * t133 / 0.2e1) * t133 + (t122 * mrSges(5,2) - t106 * mrSges(5,3) + Ifges(5,1) * t125 / 0.2e1) * t125 + (t100 * mrSges(6,1) - t101 * mrSges(6,2) + Ifges(6,3) * t123 / 0.2e1) * t123 + (t96 * mrSges(7,1) - t97 * mrSges(7,2) + Ifges(7,3) * t119 / 0.2e1) * t119 + (-t126 * mrSges(4,1) + t118 * mrSges(4,3) + Ifges(4,4) * t133 + Ifges(4,2) * t132 / 0.2e1) * t132 + (-t122 * mrSges(5,1) + t107 * mrSges(5,3) + Ifges(5,4) * t125 + Ifges(5,2) * t124 / 0.2e1) * t124 + (t104 * mrSges(6,2) - t100 * mrSges(6,3) + Ifges(6,5) * t123 + Ifges(6,1) * t121 / 0.2e1) * t121 + (t98 * mrSges(7,2) - t96 * mrSges(7,3) + Ifges(7,5) * t119 + Ifges(7,1) * t111 / 0.2e1) * t111 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t158 / 0.2e1) * t158 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t158 + Ifges(3,1) * t154 / 0.2e1) * t154) * t164 + (-t139 * t154 + t140 * t158) * mrSges(3,3)) * t164 + (t139 * mrSges(3,1) - t140 * mrSges(3,2) + Ifges(3,3) * t144 / 0.2e1 + (Ifges(3,5) * t154 + Ifges(3,6) * t158) * t164) * t144 + (-t104 * mrSges(6,1) + t101 * mrSges(6,3) + Ifges(6,4) * t121 + Ifges(6,6) * t123 + Ifges(6,2) * t120 / 0.2e1) * t120 + (-t98 * mrSges(7,1) + t97 * mrSges(7,3) + Ifges(7,4) * t111 + Ifges(7,6) * t119 + Ifges(7,2) * t110 / 0.2e1) * t110 + (t117 * mrSges(4,1) + t106 * mrSges(5,1) - t118 * mrSges(4,2) - t107 * mrSges(5,2) + Ifges(4,5) * t133 + Ifges(5,5) * t125 + Ifges(4,6) * t132 + Ifges(5,6) * t124 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t138) * t138;
T  = t1;
