% Calculate kinetic energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR15_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:27
% EndTime: 2019-03-09 20:24:27
% DurationCPUTime: 0.62s
% Computational Cost: add. (1334->120), mult. (3448->184), div. (0->0), fcn. (2832->12), ass. (0->57)
t159 = pkin(3) + pkin(11);
t158 = cos(qJ(3));
t155 = cos(pkin(6)) * qJD(1);
t134 = qJD(2) + t155;
t141 = sin(qJ(3));
t137 = cos(pkin(7));
t145 = cos(qJ(2));
t136 = sin(pkin(6));
t154 = t136 * qJD(1);
t150 = t145 * t154;
t149 = t137 * t150;
t142 = sin(qJ(2));
t151 = t142 * t154;
t135 = sin(pkin(7));
t152 = t135 * t158;
t120 = -t134 * t152 + t141 * t151 - t149 * t158;
t153 = pkin(1) * t155;
t133 = t145 * t153;
t119 = pkin(2) * t134 + t133 + (-pkin(10) * t137 - pkin(9)) * t151;
t124 = (-pkin(10) * t135 * t142 - pkin(2) * t145 - pkin(1)) * t154;
t110 = -t119 * t135 + t137 * t124;
t156 = t137 * t141;
t157 = t135 * t141;
t121 = t134 * t157 + (t142 * t158 + t145 * t156) * t154;
t148 = -qJ(4) * t121 + t110;
t101 = t120 * t159 + t148;
t140 = sin(qJ(5));
t144 = cos(qJ(5));
t126 = t134 * t137 - t135 * t150 + qJD(3);
t128 = pkin(9) * t150 + t142 * t153;
t117 = (t134 * t135 + t149) * pkin(10) + t128;
t106 = t119 * t137 * t158 - t141 * t117 + t124 * t152;
t147 = qJD(4) - t106;
t99 = t121 * pkin(4) - t126 * t159 + t147;
t96 = t144 * t101 + t140 * t99;
t107 = t158 * t117 + t119 * t156 + t124 * t157;
t105 = -t126 * qJ(4) - t107;
t95 = -t101 * t140 + t144 * t99;
t112 = t120 * t144 - t126 * t140;
t102 = -pkin(4) * t120 - t105;
t146 = qJD(1) ^ 2;
t143 = cos(qJ(6));
t139 = sin(qJ(6));
t127 = -pkin(9) * t151 + t133;
t118 = qJD(5) + t121;
t113 = t120 * t140 + t126 * t144;
t111 = qJD(6) - t112;
t109 = t113 * t143 + t118 * t139;
t108 = -t113 * t139 + t118 * t143;
t104 = -t126 * pkin(3) + t147;
t103 = pkin(3) * t120 + t148;
t97 = -pkin(5) * t112 - pkin(12) * t113 + t102;
t94 = pkin(12) * t118 + t96;
t93 = -pkin(5) * t118 - t95;
t92 = t139 * t97 + t143 * t94;
t91 = -t139 * t94 + t143 * t97;
t1 = m(5) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t146 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t136 ^ 2 * t146 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t110 ^ 2) / 0.2e1 + (t95 * mrSges(6,1) - t96 * mrSges(6,2) + Ifges(6,3) * t118 / 0.2e1) * t118 + (t91 * mrSges(7,1) - t92 * mrSges(7,2) + Ifges(7,3) * t111 / 0.2e1) * t111 + (t102 * mrSges(6,2) - t95 * mrSges(6,3) + Ifges(6,5) * t118 + Ifges(6,1) * t113 / 0.2e1) * t113 + (t93 * mrSges(7,2) - t91 * mrSges(7,3) + Ifges(7,5) * t111 + Ifges(7,1) * t109 / 0.2e1) * t109 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t145 / 0.2e1) * t145 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t145 + Ifges(3,1) * t142 / 0.2e1) * t142) * t154 + (-t127 * t142 + t128 * t145) * mrSges(3,3)) * t154 + (t127 * mrSges(3,1) - t128 * mrSges(3,2) + Ifges(3,3) * t134 / 0.2e1 + (Ifges(3,5) * t142 + Ifges(3,6) * t145) * t154) * t134 + (-t102 * mrSges(6,1) + t96 * mrSges(6,3) + Ifges(6,4) * t113 + Ifges(6,6) * t118 + Ifges(6,2) * t112 / 0.2e1) * t112 + (-t93 * mrSges(7,1) + t92 * mrSges(7,3) + Ifges(7,4) * t109 + Ifges(7,6) * t111 + Ifges(7,2) * t108 / 0.2e1) * t108 + (t106 * mrSges(4,1) - t107 * mrSges(4,2) + t104 * mrSges(5,2) - t105 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t126) * t126 + (t104 * mrSges(5,1) + t110 * mrSges(4,2) - t106 * mrSges(4,3) - t103 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t121 + (-Ifges(5,4) + Ifges(4,5)) * t126) * t121 + (t110 * mrSges(4,1) + t105 * mrSges(5,1) - t103 * mrSges(5,2) - t107 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t120 + (Ifges(5,5) - Ifges(4,6)) * t126 + (-Ifges(4,4) - Ifges(5,6)) * t121) * t120;
T  = t1;
