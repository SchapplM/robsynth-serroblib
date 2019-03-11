% Calculate kinetic energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:45
% EndTime: 2019-03-08 22:22:46
% DurationCPUTime: 0.61s
% Computational Cost: add. (796->101), mult. (1883->164), div. (0->0), fcn. (1519->14), ass. (0->51)
t152 = cos(qJ(2));
t141 = sin(pkin(6));
t160 = qJD(1) * t141;
t133 = qJD(2) * pkin(2) + t152 * t160;
t140 = sin(pkin(7));
t143 = cos(pkin(7));
t144 = cos(pkin(6));
t159 = qJD(1) * t144;
t162 = t133 * t143 + t140 * t159;
t148 = sin(qJ(2));
t158 = t140 * qJD(2);
t132 = pkin(9) * t158 + t148 * t160;
t147 = sin(qJ(3));
t151 = cos(qJ(3));
t121 = -t147 * t132 + t162 * t151;
t122 = t151 * t132 + t162 * t147;
t137 = qJD(2) * t143 + qJD(3);
t117 = qJ(4) * t137 + t122;
t136 = t143 * t159;
t125 = t136 + (-t133 + (-pkin(3) * t151 - qJ(4) * t147) * qJD(2)) * t140;
t139 = sin(pkin(13));
t142 = cos(pkin(13));
t110 = -t117 * t139 + t142 * t125;
t155 = t147 * t158;
t128 = t137 * t139 + t142 * t155;
t156 = t151 * t158;
t107 = -pkin(4) * t156 - pkin(10) * t128 + t110;
t111 = t142 * t117 + t139 * t125;
t127 = t137 * t142 - t139 * t155;
t109 = pkin(10) * t127 + t111;
t146 = sin(qJ(5));
t150 = cos(qJ(5));
t104 = t146 * t107 + t150 * t109;
t103 = t107 * t150 - t109 * t146;
t119 = t127 * t150 - t128 * t146;
t116 = -pkin(3) * t137 + qJD(4) - t121;
t112 = -pkin(4) * t127 + t116;
t149 = cos(qJ(6));
t145 = sin(qJ(6));
t135 = qJD(5) - t156;
t126 = -t133 * t140 + t136;
t120 = t127 * t146 + t128 * t150;
t118 = qJD(6) - t119;
t114 = t120 * t149 + t135 * t145;
t113 = -t120 * t145 + t135 * t149;
t105 = -pkin(5) * t119 - pkin(11) * t120 + t112;
t102 = pkin(11) * t135 + t104;
t101 = -pkin(5) * t135 - t103;
t100 = t102 * t149 + t105 * t145;
t99 = -t102 * t145 + t105 * t149;
t1 = m(5) * (t110 ^ 2 + t111 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t121 ^ 2 + t122 ^ 2 + t126 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + (t121 * mrSges(4,1) - t122 * mrSges(4,2) + Ifges(4,3) * t137 / 0.2e1) * t137 + (t103 * mrSges(6,1) - t104 * mrSges(6,2) + Ifges(6,3) * t135 / 0.2e1) * t135 + (t116 * mrSges(5,2) - t110 * mrSges(5,3) + Ifges(5,1) * t128 / 0.2e1) * t128 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t118 / 0.2e1) * t118 + (t112 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,5) * t135 + Ifges(6,1) * t120 / 0.2e1) * t120 + (t101 * mrSges(7,2) - t99 * mrSges(7,3) + Ifges(7,5) * t118 + Ifges(7,1) * t114 / 0.2e1) * t114 + (m(3) * (t144 ^ 2 + (t148 ^ 2 + t152 ^ 2) * t141 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t152 - mrSges(3,2) * t148) * t160 + ((t126 * mrSges(4,2) - t121 * mrSges(4,3) + Ifges(4,5) * t137 + Ifges(4,1) * t155 / 0.2e1) * t147 + (-t126 * mrSges(4,1) - t110 * mrSges(5,1) + t111 * mrSges(5,2) + t122 * mrSges(4,3) - Ifges(5,5) * t128 + Ifges(4,6) * t137 + (Ifges(4,4) * t147 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t151) * t158) * t151) * t140) * qJD(2) + (-Ifges(5,6) * t156 - t116 * mrSges(5,1) + t111 * mrSges(5,3) + Ifges(5,4) * t128 + Ifges(5,2) * t127 / 0.2e1) * t127 + (-t112 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t120 + Ifges(6,6) * t135 + Ifges(6,2) * t119 / 0.2e1) * t119 + (-t101 * mrSges(7,1) + t100 * mrSges(7,3) + Ifges(7,4) * t114 + Ifges(7,6) * t118 + Ifges(7,2) * t113 / 0.2e1) * t113;
T  = t1;
