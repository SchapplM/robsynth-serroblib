% Calculate kinetic energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR15_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:06
% EndTime: 2019-03-10 00:34:07
% DurationCPUTime: 0.68s
% Computational Cost: add. (1456->120), mult. (3744->184), div. (0->0), fcn. (3100->12), ass. (0->57)
t143 = sin(qJ(2));
t146 = cos(qJ(2));
t137 = sin(pkin(6));
t156 = t137 * qJD(1);
t152 = t146 * t156;
t157 = cos(pkin(6)) * qJD(1);
t155 = pkin(1) * t157;
t131 = pkin(9) * t152 + t143 * t155;
t138 = cos(pkin(7));
t135 = qJD(2) + t157;
t136 = sin(pkin(7));
t160 = t135 * t136;
t120 = (t138 * t152 + t160) * pkin(10) + t131;
t128 = (-pkin(10) * t136 * t143 - pkin(2) * t146 - pkin(1)) * t156;
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t134 = t146 * t155;
t153 = t143 * t156;
t122 = pkin(2) * t135 + t134 + (-pkin(10) * t138 - pkin(9)) * t153;
t161 = t122 * t138;
t109 = -t142 * t120 + (t128 * t136 + t161) * t145;
t158 = t138 * t146;
t123 = (-t142 * t143 + t145 * t158) * t156 + t145 * t160;
t163 = pkin(4) + pkin(12);
t162 = cos(qJ(4));
t159 = t136 * t142;
t113 = -t122 * t136 + t138 * t128;
t124 = t135 * t159 + (t142 * t158 + t143 * t145) * t156;
t104 = -pkin(3) * t123 - pkin(11) * t124 + t113;
t110 = t145 * t120 + t128 * t159 + t142 * t161;
t129 = t135 * t138 - t136 * t152 + qJD(3);
t108 = pkin(11) * t129 + t110;
t141 = sin(qJ(4));
t101 = t141 * t104 + t162 * t108;
t121 = qJD(4) - t123;
t99 = -qJ(5) * t121 - t101;
t100 = t162 * t104 - t141 * t108;
t150 = qJD(5) - t100;
t116 = t162 * t124 + t141 * t129;
t107 = -pkin(3) * t129 - t109;
t148 = -qJ(5) * t116 + t107;
t147 = qJD(1) ^ 2;
t144 = cos(qJ(6));
t140 = sin(qJ(6));
t130 = -pkin(9) * t153 + t134;
t115 = t124 * t141 - t162 * t129;
t114 = qJD(6) + t116;
t112 = t115 * t140 + t121 * t144;
t111 = t115 * t144 - t121 * t140;
t102 = pkin(4) * t115 + t148;
t98 = -t121 * pkin(4) + t150;
t97 = t163 * t115 + t148;
t96 = -pkin(5) * t115 - t99;
t95 = t116 * pkin(5) - t163 * t121 + t150;
t94 = t140 * t95 + t144 * t97;
t93 = -t140 * t97 + t144 * t95;
t1 = m(5) * (t100 ^ 2 + t101 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t93 ^ 2 + t94 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t137 ^ 2 * t147 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + t147 * Ifges(2,3) / 0.2e1 + (t109 * mrSges(4,1) - t110 * mrSges(4,2) + Ifges(4,3) * t129 / 0.2e1) * t129 + (t93 * mrSges(7,1) - t94 * mrSges(7,2) + Ifges(7,3) * t114 / 0.2e1) * t114 + (t113 * mrSges(4,2) - t109 * mrSges(4,3) + Ifges(4,5) * t129 + Ifges(4,1) * t124 / 0.2e1) * t124 + (t96 * mrSges(7,2) - t93 * mrSges(7,3) + Ifges(7,5) * t114 + Ifges(7,1) * t112 / 0.2e1) * t112 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t146 / 0.2e1) * t146 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t146 + Ifges(3,1) * t143 / 0.2e1) * t143) * t156 + (-t130 * t143 + t131 * t146) * mrSges(3,3)) * t156 + (t130 * mrSges(3,1) - t131 * mrSges(3,2) + Ifges(3,3) * t135 / 0.2e1 + (Ifges(3,5) * t143 + Ifges(3,6) * t146) * t156) * t135 + (-t113 * mrSges(4,1) + t110 * mrSges(4,3) + Ifges(4,4) * t124 + Ifges(4,6) * t129 + Ifges(4,2) * t123 / 0.2e1) * t123 + (-t96 * mrSges(7,1) + t94 * mrSges(7,3) + Ifges(7,4) * t112 + Ifges(7,6) * t114 + Ifges(7,2) * t111 / 0.2e1) * t111 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + t98 * mrSges(6,2) - t99 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t121) * t121 + (t98 * mrSges(6,1) + t107 * mrSges(5,2) - t100 * mrSges(5,3) - t102 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t116 + (-Ifges(6,4) + Ifges(5,5)) * t121) * t116 + (t107 * mrSges(5,1) + t99 * mrSges(6,1) - t102 * mrSges(6,2) - t101 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t115 + (Ifges(6,5) - Ifges(5,6)) * t121 + (-Ifges(5,4) - Ifges(6,6)) * t116) * t115;
T  = t1;
