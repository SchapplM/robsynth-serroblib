% Calculate kinetic energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:05
% EndTime: 2019-03-09 19:38:05
% DurationCPUTime: 0.66s
% Computational Cost: add. (1464->114), mult. (3258->181), div. (0->0), fcn. (2636->12), ass. (0->52)
t140 = sin(qJ(5));
t144 = cos(qJ(5));
t146 = cos(qJ(2));
t142 = sin(qJ(2));
t136 = sin(pkin(6));
t152 = qJD(1) * t136;
t149 = t142 * t152;
t151 = cos(pkin(6)) * qJD(1);
t150 = pkin(1) * t151;
t127 = -pkin(8) * t149 + t146 * t150;
t134 = qJD(2) + t151;
t120 = -pkin(2) * t134 - t127;
t141 = sin(qJ(3));
t145 = cos(qJ(3));
t125 = -t145 * t134 + t141 * t149;
t126 = t134 * t141 + t145 * t149;
t109 = pkin(3) * t125 - qJ(4) * t126 + t120;
t148 = t146 * t152;
t128 = pkin(8) * t148 + t142 * t150;
t121 = pkin(9) * t134 + t128;
t123 = (-pkin(2) * t146 - pkin(9) * t142 - pkin(1)) * t152;
t114 = t145 * t121 + t141 * t123;
t130 = qJD(3) - t148;
t112 = qJ(4) * t130 + t114;
t135 = sin(pkin(12));
t137 = cos(pkin(12));
t102 = t137 * t109 - t112 * t135;
t116 = t126 * t137 + t130 * t135;
t96 = pkin(4) * t125 - pkin(10) * t116 + t102;
t103 = t135 * t109 + t137 * t112;
t115 = -t126 * t135 + t130 * t137;
t99 = pkin(10) * t115 + t103;
t93 = t140 * t96 + t144 * t99;
t92 = -t140 * t99 + t144 * t96;
t113 = -t141 * t121 + t123 * t145;
t124 = qJD(5) + t125;
t111 = -pkin(3) * t130 + qJD(4) - t113;
t104 = -pkin(4) * t115 + t111;
t147 = qJD(1) ^ 2;
t143 = cos(qJ(6));
t139 = sin(qJ(6));
t122 = qJD(6) + t124;
t106 = t115 * t140 + t116 * t144;
t105 = t115 * t144 - t116 * t140;
t101 = t105 * t139 + t106 * t143;
t100 = t105 * t143 - t106 * t139;
t97 = -pkin(5) * t105 + t104;
t91 = pkin(11) * t105 + t93;
t90 = pkin(5) * t124 - pkin(11) * t106 + t92;
t89 = t139 * t90 + t143 * t91;
t88 = -t139 * t91 + t143 * t90;
t1 = m(4) * (t113 ^ 2 + t114 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t88 ^ 2 + t89 ^ 2 + t97 ^ 2) / 0.2e1 + t147 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t136 ^ 2 * t147 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + (t113 * mrSges(4,1) - t114 * mrSges(4,2) + Ifges(4,3) * t130 / 0.2e1) * t130 + (t92 * mrSges(6,1) - t93 * mrSges(6,2) + Ifges(6,3) * t124 / 0.2e1) * t124 + (t88 * mrSges(7,1) - t89 * mrSges(7,2) + Ifges(7,3) * t122 / 0.2e1) * t122 + (t111 * mrSges(5,2) - t102 * mrSges(5,3) + Ifges(5,1) * t116 / 0.2e1) * t116 + (t120 * mrSges(4,2) - t113 * mrSges(4,3) + Ifges(4,5) * t130 + Ifges(4,1) * t126 / 0.2e1) * t126 + (-t111 * mrSges(5,1) + t103 * mrSges(5,3) + Ifges(5,4) * t116 + Ifges(5,2) * t115 / 0.2e1) * t115 + (t104 * mrSges(6,2) - t92 * mrSges(6,3) + Ifges(6,5) * t124 + Ifges(6,1) * t106 / 0.2e1) * t106 + (t97 * mrSges(7,2) - t88 * mrSges(7,3) + Ifges(7,5) * t122 + Ifges(7,1) * t101 / 0.2e1) * t101 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t146 / 0.2e1) * t146 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t146 + Ifges(3,1) * t142 / 0.2e1) * t142) * t152 + (-t127 * t142 + t128 * t146) * mrSges(3,3)) * t152 + (t127 * mrSges(3,1) - t128 * mrSges(3,2) + Ifges(3,3) * t134 / 0.2e1 + (Ifges(3,5) * t142 + Ifges(3,6) * t146) * t152) * t134 + (-t104 * mrSges(6,1) + t93 * mrSges(6,3) + Ifges(6,4) * t106 + Ifges(6,6) * t124 + Ifges(6,2) * t105 / 0.2e1) * t105 + (-t97 * mrSges(7,1) + t89 * mrSges(7,3) + Ifges(7,4) * t101 + Ifges(7,6) * t122 + Ifges(7,2) * t100 / 0.2e1) * t100 + (t120 * mrSges(4,1) + t102 * mrSges(5,1) - t103 * mrSges(5,2) - t114 * mrSges(4,3) - Ifges(4,4) * t126 + Ifges(5,5) * t116 - Ifges(4,6) * t130 + Ifges(5,6) * t115 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t125) * t125;
T  = t1;
