% Calculate kinetic energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:36
% EndTime: 2019-03-10 04:11:37
% DurationCPUTime: 0.72s
% Computational Cost: add. (1440->114), mult. (3138->183), div. (0->0), fcn. (2538->12), ass. (0->53)
t153 = cos(pkin(6)) * qJD(1);
t136 = qJD(2) + t153;
t142 = sin(qJ(3));
t147 = cos(qJ(3));
t143 = sin(qJ(2));
t137 = sin(pkin(6));
t154 = qJD(1) * t137;
t151 = t143 * t154;
t127 = t136 * t147 - t142 * t151;
t128 = t136 * t142 + t147 * t151;
t141 = sin(qJ(4));
t146 = cos(qJ(4));
t117 = t127 * t146 - t141 * t128;
t118 = t127 * t141 + t128 * t146;
t148 = cos(qJ(2));
t152 = pkin(1) * t153;
t129 = -pkin(8) * t151 + t148 * t152;
t123 = -pkin(2) * t136 - t129;
t119 = -pkin(3) * t127 + t123;
t103 = -pkin(4) * t117 - pkin(11) * t118 + t119;
t140 = sin(qJ(5));
t145 = cos(qJ(5));
t150 = t148 * t154;
t130 = pkin(8) * t150 + t143 * t152;
t124 = pkin(9) * t136 + t130;
t126 = (-pkin(2) * t148 - pkin(9) * t143 - pkin(1)) * t154;
t113 = -t124 * t142 + t147 * t126;
t132 = qJD(3) - t150;
t107 = pkin(3) * t132 - pkin(10) * t128 + t113;
t114 = t147 * t124 + t142 * t126;
t110 = pkin(10) * t127 + t114;
t100 = t141 * t107 + t146 * t110;
t131 = qJD(4) + t132;
t98 = pkin(11) * t131 + t100;
t94 = t140 * t103 + t145 * t98;
t93 = t145 * t103 - t140 * t98;
t99 = t107 * t146 - t141 * t110;
t116 = qJD(5) - t117;
t97 = -pkin(4) * t131 - t99;
t149 = qJD(1) ^ 2;
t144 = cos(qJ(6));
t139 = sin(qJ(6));
t115 = qJD(6) + t116;
t112 = t118 * t145 + t131 * t140;
t111 = -t118 * t140 + t131 * t145;
t105 = t111 * t139 + t112 * t144;
t104 = t111 * t144 - t112 * t139;
t95 = -pkin(5) * t111 + t97;
t92 = pkin(12) * t111 + t94;
t91 = pkin(5) * t116 - pkin(12) * t112 + t93;
t90 = t139 * t91 + t144 * t92;
t89 = -t139 * t92 + t144 * t91;
t1 = m(7) * (t89 ^ 2 + t90 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t119 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t123 ^ 2) / 0.2e1 + t149 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t137 ^ 2 * t149 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + (t113 * mrSges(4,1) - t114 * mrSges(4,2) + Ifges(4,3) * t132 / 0.2e1) * t132 + (t99 * mrSges(5,1) - t100 * mrSges(5,2) + Ifges(5,3) * t131 / 0.2e1) * t131 + (t93 * mrSges(6,1) - t94 * mrSges(6,2) + Ifges(6,3) * t116 / 0.2e1) * t116 + (t89 * mrSges(7,1) - t90 * mrSges(7,2) + Ifges(7,3) * t115 / 0.2e1) * t115 + (t123 * mrSges(4,2) - t113 * mrSges(4,3) + Ifges(4,5) * t132 + Ifges(4,1) * t128 / 0.2e1) * t128 + (t119 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,5) * t131 + Ifges(5,1) * t118 / 0.2e1) * t118 + (t97 * mrSges(6,2) - t93 * mrSges(6,3) + Ifges(6,5) * t116 + Ifges(6,1) * t112 / 0.2e1) * t112 + (t95 * mrSges(7,2) - t89 * mrSges(7,3) + Ifges(7,5) * t115 + Ifges(7,1) * t105 / 0.2e1) * t105 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t148 / 0.2e1) * t148 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t148 + Ifges(3,1) * t143 / 0.2e1) * t143) * t154 + (-t129 * t143 + t130 * t148) * mrSges(3,3)) * t154 + (t129 * mrSges(3,1) - t130 * mrSges(3,2) + Ifges(3,3) * t136 / 0.2e1 + (Ifges(3,5) * t143 + Ifges(3,6) * t148) * t154) * t136 + (-t123 * mrSges(4,1) + t114 * mrSges(4,3) + Ifges(4,4) * t128 + Ifges(4,6) * t132 + Ifges(4,2) * t127 / 0.2e1) * t127 + (-t119 * mrSges(5,1) + t100 * mrSges(5,3) + Ifges(5,4) * t118 + Ifges(5,6) * t131 + Ifges(5,2) * t117 / 0.2e1) * t117 + (-t97 * mrSges(6,1) + t94 * mrSges(6,3) + Ifges(6,4) * t112 + Ifges(6,6) * t116 + Ifges(6,2) * t111 / 0.2e1) * t111 + (-t95 * mrSges(7,1) + t90 * mrSges(7,3) + Ifges(7,4) * t105 + Ifges(7,6) * t115 + Ifges(7,2) * t104 / 0.2e1) * t104;
T  = t1;
