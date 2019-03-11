% Calculate kinetic energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:45
% EndTime: 2019-03-09 00:54:46
% DurationCPUTime: 0.57s
% Computational Cost: add. (828->101), mult. (1883->169), div. (0->0), fcn. (1519->14), ass. (0->51)
t151 = cos(qJ(2));
t139 = sin(pkin(6));
t158 = qJD(1) * t139;
t131 = qJD(2) * pkin(2) + t151 * t158;
t138 = sin(pkin(7));
t140 = cos(pkin(7));
t141 = cos(pkin(6));
t157 = qJD(1) * t141;
t160 = t131 * t140 + t138 * t157;
t146 = sin(qJ(2));
t156 = qJD(2) * t138;
t130 = pkin(9) * t156 + t146 * t158;
t145 = sin(qJ(3));
t150 = cos(qJ(3));
t119 = -t130 * t145 + t150 * t160;
t120 = t150 * t130 + t145 * t160;
t136 = qJD(2) * t140 + qJD(3);
t115 = pkin(10) * t136 + t120;
t135 = t140 * t157;
t123 = t135 + (-t131 + (-pkin(3) * t150 - pkin(10) * t145) * qJD(2)) * t138;
t144 = sin(qJ(4));
t149 = cos(qJ(4));
t108 = -t115 * t144 + t123 * t149;
t154 = t145 * t156;
t126 = t136 * t144 + t149 * t154;
t134 = -t150 * t156 + qJD(4);
t105 = pkin(4) * t134 - pkin(11) * t126 + t108;
t109 = t115 * t149 + t123 * t144;
t125 = t136 * t149 - t144 * t154;
t107 = pkin(11) * t125 + t109;
t143 = sin(qJ(5));
t148 = cos(qJ(5));
t102 = t105 * t143 + t107 * t148;
t101 = t105 * t148 - t107 * t143;
t117 = t125 * t148 - t126 * t143;
t114 = -pkin(3) * t136 - t119;
t110 = -pkin(4) * t125 + t114;
t147 = cos(qJ(6));
t142 = sin(qJ(6));
t132 = qJD(5) + t134;
t124 = -t131 * t138 + t135;
t118 = t125 * t143 + t126 * t148;
t116 = qJD(6) - t117;
t112 = t118 * t147 + t132 * t142;
t111 = -t118 * t142 + t132 * t147;
t103 = -pkin(5) * t117 - pkin(12) * t118 + t110;
t100 = pkin(12) * t132 + t102;
t99 = -pkin(5) * t132 - t101;
t98 = t100 * t147 + t103 * t142;
t97 = -t100 * t142 + t103 * t147;
t1 = m(4) * (t119 ^ 2 + t120 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t114 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + (t119 * mrSges(4,1) - t120 * mrSges(4,2) + Ifges(4,3) * t136 / 0.2e1) * t136 + (t108 * mrSges(5,1) - t109 * mrSges(5,2) + Ifges(5,3) * t134 / 0.2e1) * t134 + (t101 * mrSges(6,1) - t102 * mrSges(6,2) + Ifges(6,3) * t132 / 0.2e1) * t132 + (t97 * mrSges(7,1) - t98 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t114 * mrSges(5,2) - t108 * mrSges(5,3) + Ifges(5,5) * t134 + Ifges(5,1) * t126 / 0.2e1) * t126 + (t110 * mrSges(6,2) - t101 * mrSges(6,3) + Ifges(6,5) * t132 + Ifges(6,1) * t118 / 0.2e1) * t118 + (t99 * mrSges(7,2) - t97 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t112 / 0.2e1) * t112 + (m(2) / 0.2e1 + m(3) * (t141 ^ 2 + (t146 ^ 2 + t151 ^ 2) * t139 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t114 * mrSges(5,1) + t109 * mrSges(5,3) + Ifges(5,4) * t126 + Ifges(5,6) * t134 + Ifges(5,2) * t125 / 0.2e1) * t125 + (-t110 * mrSges(6,1) + t102 * mrSges(6,3) + Ifges(6,4) * t118 + Ifges(6,6) * t132 + Ifges(6,2) * t117 / 0.2e1) * t117 + (-t99 * mrSges(7,1) + t98 * mrSges(7,3) + Ifges(7,4) * t112 + Ifges(7,6) * t116 + Ifges(7,2) * t111 / 0.2e1) * t111 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t151 - mrSges(3,2) * t146) * t158 + (t124 * (-mrSges(4,1) * t150 + mrSges(4,2) * t145) + (Ifges(4,2) * t150 ^ 2 / 0.2e1 + (Ifges(4,4) * t150 + Ifges(4,1) * t145 / 0.2e1) * t145) * t156 + (-t119 * t145 + t120 * t150) * mrSges(4,3) + t136 * (Ifges(4,5) * t145 + Ifges(4,6) * t150)) * t138) * qJD(2);
T  = t1;
