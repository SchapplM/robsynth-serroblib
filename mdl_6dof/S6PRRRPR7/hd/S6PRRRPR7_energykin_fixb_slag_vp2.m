% Calculate kinetic energy for
% S6PRRRPR7
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:38:43
% EndTime: 2019-03-08 23:38:44
% DurationCPUTime: 0.54s
% Computational Cost: add. (810->101), mult. (1857->167), div. (0->0), fcn. (1495->14), ass. (0->50)
t152 = cos(qJ(2));
t142 = sin(pkin(6));
t159 = qJD(1) * t142;
t133 = qJD(2) * pkin(2) + t152 * t159;
t141 = sin(pkin(7));
t144 = cos(pkin(7));
t145 = cos(pkin(6));
t158 = qJD(1) * t145;
t162 = t133 * t144 + t141 * t158;
t149 = sin(qJ(2));
t157 = qJD(2) * t141;
t132 = pkin(9) * t157 + t149 * t159;
t148 = sin(qJ(3));
t151 = cos(qJ(3));
t119 = -t148 * t132 + t162 * t151;
t161 = cos(qJ(4));
t120 = t151 * t132 + t162 * t148;
t138 = qJD(2) * t144 + qJD(3);
t118 = pkin(10) * t138 + t120;
t137 = t144 * t158;
t122 = t137 + (-t133 + (-pkin(3) * t151 - pkin(10) * t148) * qJD(2)) * t141;
t147 = sin(qJ(4));
t111 = t161 * t118 + t147 * t122;
t136 = -t151 * t157 + qJD(4);
t107 = qJ(5) * t136 + t111;
t117 = -pkin(3) * t138 - t119;
t155 = t148 * t157;
t127 = -t161 * t138 + t147 * t155;
t128 = t147 * t138 + t161 * t155;
t112 = pkin(4) * t127 - qJ(5) * t128 + t117;
t140 = sin(pkin(13));
t143 = cos(pkin(13));
t103 = t143 * t107 + t140 * t112;
t102 = -t107 * t140 + t143 * t112;
t110 = -t147 * t118 + t161 * t122;
t106 = -t136 * pkin(4) + qJD(5) - t110;
t150 = cos(qJ(6));
t146 = sin(qJ(6));
t126 = qJD(6) + t127;
t125 = -t133 * t141 + t137;
t124 = t128 * t143 + t136 * t140;
t123 = -t128 * t140 + t136 * t143;
t114 = t123 * t146 + t124 * t150;
t113 = t123 * t150 - t124 * t146;
t104 = -t123 * pkin(5) + t106;
t101 = pkin(11) * t123 + t103;
t100 = pkin(5) * t127 - pkin(11) * t124 + t102;
t99 = t100 * t146 + t101 * t150;
t98 = t100 * t150 - t101 * t146;
t1 = m(4) * (t119 ^ 2 + t120 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t106 ^ 2) / 0.2e1 + (t119 * mrSges(4,1) - t120 * mrSges(4,2) + Ifges(4,3) * t138 / 0.2e1) * t138 + (t110 * mrSges(5,1) - t111 * mrSges(5,2) + Ifges(5,3) * t136 / 0.2e1) * t136 + (t98 * mrSges(7,1) - t99 * mrSges(7,2) + Ifges(7,3) * t126 / 0.2e1) * t126 + (t106 * mrSges(6,2) - t102 * mrSges(6,3) + Ifges(6,1) * t124 / 0.2e1) * t124 + (t117 * mrSges(5,2) - t110 * mrSges(5,3) + Ifges(5,5) * t136 + Ifges(5,1) * t128 / 0.2e1) * t128 + (-t106 * mrSges(6,1) + t103 * mrSges(6,3) + Ifges(6,4) * t124 + Ifges(6,2) * t123 / 0.2e1) * t123 + (t104 * mrSges(7,2) - t98 * mrSges(7,3) + Ifges(7,5) * t126 + Ifges(7,1) * t114 / 0.2e1) * t114 + (m(3) * (t145 ^ 2 + (t149 ^ 2 + t152 ^ 2) * t142 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t104 * mrSges(7,1) + t99 * mrSges(7,3) + Ifges(7,4) * t114 + Ifges(7,6) * t126 + Ifges(7,2) * t113 / 0.2e1) * t113 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t152 - mrSges(3,2) * t149) * t159 + (t125 * (-mrSges(4,1) * t151 + mrSges(4,2) * t148) + (Ifges(4,2) * t151 ^ 2 / 0.2e1 + (Ifges(4,4) * t151 + Ifges(4,1) * t148 / 0.2e1) * t148) * t157 + (-t119 * t148 + t120 * t151) * mrSges(4,3) + t138 * (Ifges(4,5) * t148 + Ifges(4,6) * t151)) * t141) * qJD(2) + (t117 * mrSges(5,1) + t102 * mrSges(6,1) - t103 * mrSges(6,2) - t111 * mrSges(5,3) - Ifges(5,4) * t128 + Ifges(6,5) * t124 - Ifges(5,6) * t136 + Ifges(6,6) * t123 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t127) * t127;
T  = t1;
