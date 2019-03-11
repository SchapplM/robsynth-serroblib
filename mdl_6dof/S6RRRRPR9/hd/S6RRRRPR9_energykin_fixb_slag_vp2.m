% Calculate kinetic energy for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:42
% EndTime: 2019-03-09 22:47:42
% DurationCPUTime: 0.66s
% Computational Cost: add. (1428->114), mult. (3138->181), div. (0->0), fcn. (2538->12), ass. (0->52)
t152 = cos(qJ(4));
t150 = cos(pkin(6)) * qJD(1);
t134 = qJD(2) + t150;
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t142 = sin(qJ(2));
t136 = sin(pkin(6));
t151 = qJD(1) * t136;
t148 = t142 * t151;
t125 = t134 * t144 - t141 * t148;
t126 = t134 * t141 + t144 * t148;
t140 = sin(qJ(4));
t115 = -t152 * t125 + t126 * t140;
t116 = t140 * t125 + t152 * t126;
t145 = cos(qJ(2));
t149 = pkin(1) * t150;
t127 = -pkin(8) * t148 + t145 * t149;
t121 = -pkin(2) * t134 - t127;
t117 = -pkin(3) * t125 + t121;
t102 = pkin(4) * t115 - qJ(5) * t116 + t117;
t135 = sin(pkin(12));
t137 = cos(pkin(12));
t147 = t145 * t151;
t130 = qJD(3) - t147;
t129 = qJD(4) + t130;
t128 = pkin(8) * t147 + t142 * t149;
t122 = pkin(9) * t134 + t128;
t124 = (-pkin(2) * t145 - pkin(9) * t142 - pkin(1)) * t151;
t112 = -t122 * t141 + t144 * t124;
t106 = pkin(3) * t130 - pkin(10) * t126 + t112;
t113 = t144 * t122 + t141 * t124;
t109 = pkin(10) * t125 + t113;
t99 = t140 * t106 + t152 * t109;
t97 = qJ(5) * t129 + t99;
t93 = t135 * t102 + t137 * t97;
t92 = t137 * t102 - t135 * t97;
t98 = t152 * t106 - t140 * t109;
t96 = -t129 * pkin(4) + qJD(5) - t98;
t146 = qJD(1) ^ 2;
t143 = cos(qJ(6));
t139 = sin(qJ(6));
t114 = qJD(6) + t115;
t111 = t116 * t137 + t129 * t135;
t110 = -t116 * t135 + t129 * t137;
t104 = t110 * t139 + t111 * t143;
t103 = t110 * t143 - t111 * t139;
t94 = -t110 * pkin(5) + t96;
t91 = pkin(11) * t110 + t93;
t90 = pkin(5) * t115 - pkin(11) * t111 + t92;
t89 = t139 * t90 + t143 * t91;
t88 = -t139 * t91 + t143 * t90;
t1 = m(5) * (t117 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t121 ^ 2) / 0.2e1 + m(7) * (t88 ^ 2 + t89 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t92 ^ 2 + t93 ^ 2 + t96 ^ 2) / 0.2e1 + t146 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t136 ^ 2 * t146 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + (t112 * mrSges(4,1) - t113 * mrSges(4,2) + Ifges(4,3) * t130 / 0.2e1) * t130 + (t98 * mrSges(5,1) - t99 * mrSges(5,2) + Ifges(5,3) * t129 / 0.2e1) * t129 + (t88 * mrSges(7,1) - t89 * mrSges(7,2) + Ifges(7,3) * t114 / 0.2e1) * t114 + (t96 * mrSges(6,2) - t92 * mrSges(6,3) + Ifges(6,1) * t111 / 0.2e1) * t111 + (t121 * mrSges(4,2) - t112 * mrSges(4,3) + Ifges(4,5) * t130 + Ifges(4,1) * t126 / 0.2e1) * t126 + (t117 * mrSges(5,2) - t98 * mrSges(5,3) + Ifges(5,5) * t129 + Ifges(5,1) * t116 / 0.2e1) * t116 + (-t96 * mrSges(6,1) + t93 * mrSges(6,3) + Ifges(6,4) * t111 + Ifges(6,2) * t110 / 0.2e1) * t110 + (t94 * mrSges(7,2) - t88 * mrSges(7,3) + Ifges(7,5) * t114 + Ifges(7,1) * t104 / 0.2e1) * t104 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t145 / 0.2e1) * t145 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t145 + Ifges(3,1) * t142 / 0.2e1) * t142) * t151 + (-t127 * t142 + t128 * t145) * mrSges(3,3)) * t151 + (t127 * mrSges(3,1) - t128 * mrSges(3,2) + Ifges(3,3) * t134 / 0.2e1 + (Ifges(3,5) * t142 + Ifges(3,6) * t145) * t151) * t134 + (-t121 * mrSges(4,1) + t113 * mrSges(4,3) + Ifges(4,4) * t126 + Ifges(4,6) * t130 + Ifges(4,2) * t125 / 0.2e1) * t125 + (-t94 * mrSges(7,1) + t89 * mrSges(7,3) + Ifges(7,4) * t104 + Ifges(7,6) * t114 + Ifges(7,2) * t103 / 0.2e1) * t103 + (t117 * mrSges(5,1) + t92 * mrSges(6,1) - t93 * mrSges(6,2) - t99 * mrSges(5,3) - Ifges(5,4) * t116 + Ifges(6,5) * t111 - Ifges(5,6) * t129 + Ifges(6,6) * t110 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t115) * t115;
T  = t1;
