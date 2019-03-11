% Calculate kinetic energy for
% S6PRRRRR5
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:49
% EndTime: 2019-03-09 01:03:50
% DurationCPUTime: 0.57s
% Computational Cost: add. (822->101), mult. (1857->169), div. (0->0), fcn. (1495->14), ass. (0->51)
t155 = cos(qJ(2));
t143 = sin(pkin(6));
t162 = qJD(1) * t143;
t135 = qJD(2) * pkin(2) + t155 * t162;
t142 = sin(pkin(7));
t144 = cos(pkin(7));
t145 = cos(pkin(6));
t161 = qJD(1) * t145;
t164 = t135 * t144 + t142 * t161;
t150 = sin(qJ(2));
t160 = qJD(2) * t142;
t134 = pkin(9) * t160 + t150 * t162;
t149 = sin(qJ(3));
t154 = cos(qJ(3));
t120 = -t134 * t149 + t154 * t164;
t121 = t154 * t134 + t149 * t164;
t140 = qJD(2) * t144 + qJD(3);
t119 = pkin(10) * t140 + t121;
t139 = t144 * t161;
t123 = t139 + (-t135 + (-pkin(3) * t154 - pkin(10) * t149) * qJD(2)) * t142;
t148 = sin(qJ(4));
t153 = cos(qJ(4));
t112 = t119 * t153 + t123 * t148;
t138 = -t154 * t160 + qJD(4);
t108 = pkin(11) * t138 + t112;
t118 = -pkin(3) * t140 - t120;
t158 = t149 * t160;
t129 = t140 * t153 - t148 * t158;
t130 = t140 * t148 + t153 * t158;
t113 = -pkin(4) * t129 - pkin(11) * t130 + t118;
t147 = sin(qJ(5));
t152 = cos(qJ(5));
t104 = t108 * t152 + t113 * t147;
t103 = -t108 * t147 + t113 * t152;
t111 = -t119 * t148 + t123 * t153;
t128 = qJD(5) - t129;
t107 = -pkin(4) * t138 - t111;
t151 = cos(qJ(6));
t146 = sin(qJ(6));
t127 = -t135 * t142 + t139;
t126 = qJD(6) + t128;
t125 = t130 * t152 + t138 * t147;
t124 = -t130 * t147 + t138 * t152;
t115 = t124 * t146 + t125 * t151;
t114 = t124 * t151 - t125 * t146;
t105 = -pkin(5) * t124 + t107;
t102 = pkin(12) * t124 + t104;
t101 = pkin(5) * t128 - pkin(12) * t125 + t103;
t100 = t101 * t146 + t102 * t151;
t99 = t101 * t151 - t102 * t146;
t1 = m(4) * (t120 ^ 2 + t121 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t107 ^ 2) / 0.2e1 + (t120 * mrSges(4,1) - t121 * mrSges(4,2) + Ifges(4,3) * t140 / 0.2e1) * t140 + (t111 * mrSges(5,1) - t112 * mrSges(5,2) + Ifges(5,3) * t138 / 0.2e1) * t138 + (t103 * mrSges(6,1) - t104 * mrSges(6,2) + Ifges(6,3) * t128 / 0.2e1) * t128 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t126 / 0.2e1) * t126 + (t118 * mrSges(5,2) - t111 * mrSges(5,3) + Ifges(5,5) * t138 + Ifges(5,1) * t130 / 0.2e1) * t130 + (t107 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,5) * t128 + Ifges(6,1) * t125 / 0.2e1) * t125 + (t105 * mrSges(7,2) - t99 * mrSges(7,3) + Ifges(7,5) * t126 + Ifges(7,1) * t115 / 0.2e1) * t115 + (m(3) * (t145 ^ 2 + (t150 ^ 2 + t155 ^ 2) * t143 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t118 * mrSges(5,1) + t112 * mrSges(5,3) + Ifges(5,4) * t130 + Ifges(5,6) * t138 + Ifges(5,2) * t129 / 0.2e1) * t129 + (-t107 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t125 + Ifges(6,6) * t128 + Ifges(6,2) * t124 / 0.2e1) * t124 + (-t105 * mrSges(7,1) + t100 * mrSges(7,3) + Ifges(7,4) * t115 + Ifges(7,6) * t126 + Ifges(7,2) * t114 / 0.2e1) * t114 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t155 - mrSges(3,2) * t150) * t162 + (t127 * (-mrSges(4,1) * t154 + mrSges(4,2) * t149) + (Ifges(4,2) * t154 ^ 2 / 0.2e1 + (Ifges(4,4) * t154 + Ifges(4,1) * t149 / 0.2e1) * t149) * t160 + (-t120 * t149 + t121 * t154) * mrSges(4,3) + t140 * (Ifges(4,5) * t149 + Ifges(4,6) * t154)) * t142) * qJD(2);
T  = t1;
