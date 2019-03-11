% Calculate kinetic energy for
% S6RRRRPR12
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:15
% EndTime: 2019-03-09 23:29:16
% DurationCPUTime: 0.75s
% Computational Cost: add. (2196->123), mult. (5644->199), div. (0->0), fcn. (4776->14), ass. (0->60)
t154 = sin(qJ(2));
t158 = cos(qJ(2));
t147 = sin(pkin(6));
t167 = qJD(1) * t147;
t162 = t158 * t167;
t166 = cos(pkin(6)) * qJD(1);
t165 = pkin(1) * t166;
t140 = pkin(9) * t162 + t154 * t165;
t149 = cos(pkin(7));
t144 = qJD(2) + t166;
t146 = sin(pkin(7));
t170 = t144 * t146;
t130 = (t149 * t162 + t170) * pkin(10) + t140;
t137 = (-pkin(10) * t146 * t154 - pkin(2) * t158 - pkin(1)) * t167;
t153 = sin(qJ(3));
t157 = cos(qJ(3));
t143 = t158 * t165;
t163 = t154 * t167;
t132 = pkin(2) * t144 + t143 + (-pkin(10) * t149 - pkin(9)) * t163;
t171 = t132 * t149;
t122 = -t153 * t130 + (t137 * t146 + t171) * t157;
t168 = t149 * t158;
t133 = (-t153 * t154 + t157 * t168) * t167 + t157 * t170;
t169 = t146 * t153;
t124 = -t132 * t146 + t149 * t137;
t134 = t144 * t169 + (t153 * t168 + t154 * t157) * t167;
t115 = -pkin(3) * t133 - pkin(11) * t134 + t124;
t123 = t157 * t130 + t137 * t169 + t153 * t171;
t138 = t144 * t149 - t146 * t162 + qJD(3);
t118 = pkin(11) * t138 + t123;
t152 = sin(qJ(4));
t156 = cos(qJ(4));
t108 = t156 * t115 - t118 * t152;
t126 = t134 * t156 + t138 * t152;
t131 = qJD(4) - t133;
t105 = pkin(4) * t131 - qJ(5) * t126 + t108;
t109 = t152 * t115 + t156 * t118;
t125 = -t134 * t152 + t138 * t156;
t107 = qJ(5) * t125 + t109;
t145 = sin(pkin(13));
t148 = cos(pkin(13));
t102 = t145 * t105 + t148 * t107;
t101 = t105 * t148 - t107 * t145;
t120 = t125 * t148 - t126 * t145;
t117 = -pkin(3) * t138 - t122;
t110 = -pkin(4) * t125 + qJD(5) + t117;
t159 = qJD(1) ^ 2;
t155 = cos(qJ(6));
t151 = sin(qJ(6));
t139 = -pkin(9) * t163 + t143;
t121 = t125 * t145 + t126 * t148;
t119 = qJD(6) - t120;
t112 = t121 * t155 + t131 * t151;
t111 = -t121 * t151 + t131 * t155;
t103 = -pkin(5) * t120 - pkin(12) * t121 + t110;
t100 = pkin(12) * t131 + t102;
t99 = -pkin(5) * t131 - t101;
t98 = t100 * t155 + t103 * t151;
t97 = -t100 * t151 + t103 * t155;
t1 = m(7) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t117 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + t159 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t147 ^ 2 * t159 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + (t122 * mrSges(4,1) - t123 * mrSges(4,2) + Ifges(4,3) * t138 / 0.2e1) * t138 + (t117 * mrSges(5,2) - t108 * mrSges(5,3) + Ifges(5,1) * t126 / 0.2e1) * t126 + (t110 * mrSges(6,2) - t101 * mrSges(6,3) + Ifges(6,1) * t121 / 0.2e1) * t121 + (t97 * mrSges(7,1) - t98 * mrSges(7,2) + Ifges(7,3) * t119 / 0.2e1) * t119 + (t124 * mrSges(4,2) - t122 * mrSges(4,3) + Ifges(4,5) * t138 + Ifges(4,1) * t134 / 0.2e1) * t134 + (-t117 * mrSges(5,1) + t109 * mrSges(5,3) + Ifges(5,4) * t126 + Ifges(5,2) * t125 / 0.2e1) * t125 + (-t110 * mrSges(6,1) + t102 * mrSges(6,3) + Ifges(6,4) * t121 + Ifges(6,2) * t120 / 0.2e1) * t120 + (t99 * mrSges(7,2) - t97 * mrSges(7,3) + Ifges(7,5) * t119 + Ifges(7,1) * t112 / 0.2e1) * t112 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t158 / 0.2e1) * t158 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t158 + Ifges(3,1) * t154 / 0.2e1) * t154) * t167 + (-t139 * t154 + t140 * t158) * mrSges(3,3)) * t167 + (t139 * mrSges(3,1) - t140 * mrSges(3,2) + Ifges(3,3) * t144 / 0.2e1 + (Ifges(3,5) * t154 + Ifges(3,6) * t158) * t167) * t144 + (-t124 * mrSges(4,1) + t123 * mrSges(4,3) + Ifges(4,4) * t134 + Ifges(4,6) * t138 + Ifges(4,2) * t133 / 0.2e1) * t133 + (-t99 * mrSges(7,1) + t98 * mrSges(7,3) + Ifges(7,4) * t112 + Ifges(7,6) * t119 + Ifges(7,2) * t111 / 0.2e1) * t111 + (t108 * mrSges(5,1) + t101 * mrSges(6,1) - t109 * mrSges(5,2) - t102 * mrSges(6,2) + Ifges(5,5) * t126 + Ifges(6,5) * t121 + Ifges(5,6) * t125 + Ifges(6,6) * t120 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t131) * t131;
T  = t1;
