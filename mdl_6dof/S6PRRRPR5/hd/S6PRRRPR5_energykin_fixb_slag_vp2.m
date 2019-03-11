% Calculate kinetic energy for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:22:54
% EndTime: 2019-03-08 23:22:55
% DurationCPUTime: 0.52s
% Computational Cost: add. (810->101), mult. (1883->167), div. (0->0), fcn. (1519->14), ass. (0->50)
t149 = cos(qJ(2));
t138 = sin(pkin(6));
t156 = qJD(1) * t138;
t130 = qJD(2) * pkin(2) + t149 * t156;
t137 = sin(pkin(7));
t140 = cos(pkin(7));
t141 = cos(pkin(6));
t155 = qJD(1) * t141;
t158 = t130 * t140 + t137 * t155;
t145 = sin(qJ(2));
t154 = qJD(2) * t137;
t129 = pkin(9) * t154 + t145 * t156;
t144 = sin(qJ(3));
t148 = cos(qJ(3));
t118 = -t144 * t129 + t158 * t148;
t119 = t148 * t129 + t158 * t144;
t134 = qJD(2) * t140 + qJD(3);
t114 = pkin(10) * t134 + t119;
t133 = t140 * t155;
t122 = t133 + (-t130 + (-pkin(3) * t148 - pkin(10) * t144) * qJD(2)) * t137;
t143 = sin(qJ(4));
t147 = cos(qJ(4));
t107 = -t114 * t143 + t147 * t122;
t152 = t144 * t154;
t125 = t134 * t143 + t147 * t152;
t132 = -t148 * t154 + qJD(4);
t104 = pkin(4) * t132 - qJ(5) * t125 + t107;
t108 = t147 * t114 + t143 * t122;
t124 = t134 * t147 - t143 * t152;
t106 = qJ(5) * t124 + t108;
t136 = sin(pkin(13));
t139 = cos(pkin(13));
t101 = t136 * t104 + t139 * t106;
t100 = t104 * t139 - t106 * t136;
t116 = t124 * t139 - t125 * t136;
t113 = -pkin(3) * t134 - t118;
t109 = -pkin(4) * t124 + qJD(5) + t113;
t146 = cos(qJ(6));
t142 = sin(qJ(6));
t123 = -t130 * t137 + t133;
t117 = t124 * t136 + t125 * t139;
t115 = qJD(6) - t116;
t111 = t117 * t146 + t132 * t142;
t110 = -t117 * t142 + t132 * t146;
t102 = -pkin(5) * t116 - pkin(11) * t117 + t109;
t99 = pkin(11) * t132 + t101;
t98 = -pkin(5) * t132 - t100;
t97 = t102 * t142 + t146 * t99;
t96 = t102 * t146 - t142 * t99;
t1 = m(4) * (t118 ^ 2 + t119 ^ 2 + t123 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + (t118 * mrSges(4,1) - t119 * mrSges(4,2) + Ifges(4,3) * t134 / 0.2e1) * t134 + (t113 * mrSges(5,2) - t107 * mrSges(5,3) + Ifges(5,1) * t125 / 0.2e1) * t125 + (t109 * mrSges(6,2) - t100 * mrSges(6,3) + Ifges(6,1) * t117 / 0.2e1) * t117 + (t96 * mrSges(7,1) - t97 * mrSges(7,2) + Ifges(7,3) * t115 / 0.2e1) * t115 + (-t113 * mrSges(5,1) + t108 * mrSges(5,3) + Ifges(5,4) * t125 + Ifges(5,2) * t124 / 0.2e1) * t124 + (-t109 * mrSges(6,1) + t101 * mrSges(6,3) + Ifges(6,4) * t117 + Ifges(6,2) * t116 / 0.2e1) * t116 + (t98 * mrSges(7,2) - t96 * mrSges(7,3) + Ifges(7,5) * t115 + Ifges(7,1) * t111 / 0.2e1) * t111 + (m(3) * (t141 ^ 2 + (t145 ^ 2 + t149 ^ 2) * t138 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t98 * mrSges(7,1) + t97 * mrSges(7,3) + Ifges(7,4) * t111 + Ifges(7,6) * t115 + Ifges(7,2) * t110 / 0.2e1) * t110 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t149 - mrSges(3,2) * t145) * t156 + (t123 * (-mrSges(4,1) * t148 + mrSges(4,2) * t144) + (Ifges(4,2) * t148 ^ 2 / 0.2e1 + (Ifges(4,4) * t148 + Ifges(4,1) * t144 / 0.2e1) * t144) * t154 + (-t118 * t144 + t119 * t148) * mrSges(4,3) + t134 * (Ifges(4,5) * t144 + Ifges(4,6) * t148)) * t137) * qJD(2) + (t107 * mrSges(5,1) + t100 * mrSges(6,1) - t108 * mrSges(5,2) - t101 * mrSges(6,2) + Ifges(5,5) * t125 + Ifges(6,5) * t117 + Ifges(5,6) * t124 + Ifges(6,6) * t116 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t132) * t132;
T  = t1;
