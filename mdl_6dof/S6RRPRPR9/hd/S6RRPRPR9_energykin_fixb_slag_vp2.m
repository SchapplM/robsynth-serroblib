% Calculate kinetic energy for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:19
% EndTime: 2019-03-09 10:55:20
% DurationCPUTime: 0.61s
% Computational Cost: add. (1350->114), mult. (3138->178), div. (0->0), fcn. (2538->12), ass. (0->51)
t153 = cos(qJ(4));
t146 = cos(qJ(2));
t144 = sin(qJ(2));
t138 = sin(pkin(6));
t151 = t138 * qJD(1);
t149 = t144 * t151;
t152 = cos(pkin(6)) * qJD(1);
t150 = pkin(1) * t152;
t129 = -pkin(8) * t149 + t146 * t150;
t135 = qJD(2) + t152;
t121 = -pkin(2) * t135 + qJD(3) - t129;
t137 = sin(pkin(11));
t140 = cos(pkin(11));
t127 = t135 * t140 - t137 * t149;
t117 = -pkin(3) * t127 + t121;
t128 = t135 * t137 + t140 * t149;
t143 = sin(qJ(4));
t118 = -t153 * t127 + t128 * t143;
t119 = t143 * t127 + t128 * t153;
t104 = pkin(4) * t118 - qJ(5) * t119 + t117;
t136 = sin(pkin(12));
t139 = cos(pkin(12));
t148 = t146 * t151;
t130 = pkin(8) * t148 + t144 * t150;
t124 = qJ(3) * t135 + t130;
t126 = (-pkin(2) * t146 - qJ(3) * t144 - pkin(1)) * t151;
t114 = -t124 * t137 + t140 * t126;
t108 = -pkin(3) * t148 - pkin(9) * t128 + t114;
t115 = t140 * t124 + t137 * t126;
t111 = pkin(9) * t127 + t115;
t101 = t143 * t108 + t153 * t111;
t131 = qJD(4) - t148;
t99 = qJ(5) * t131 + t101;
t95 = t136 * t104 + t139 * t99;
t94 = t139 * t104 - t136 * t99;
t100 = t108 * t153 - t143 * t111;
t98 = -t131 * pkin(4) + qJD(5) - t100;
t147 = qJD(1) ^ 2;
t145 = cos(qJ(6));
t142 = sin(qJ(6));
t116 = qJD(6) + t118;
t113 = t119 * t139 + t131 * t136;
t112 = -t119 * t136 + t131 * t139;
t106 = t112 * t142 + t113 * t145;
t105 = t112 * t145 - t113 * t142;
t96 = -t112 * pkin(5) + t98;
t93 = pkin(10) * t112 + t95;
t92 = pkin(5) * t118 - pkin(10) * t113 + t94;
t91 = t142 * t92 + t145 * t93;
t90 = -t142 * t93 + t145 * t92;
t1 = m(3) * (pkin(1) ^ 2 * t138 ^ 2 * t147 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + t147 * Ifges(2,3) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) / 0.2e1 + (t129 * mrSges(3,1) - t130 * mrSges(3,2) + Ifges(3,3) * t135 / 0.2e1) * t135 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t131 / 0.2e1) * t131 + (t121 * mrSges(4,2) - t114 * mrSges(4,3) + Ifges(4,1) * t128 / 0.2e1) * t128 + (t90 * mrSges(7,1) - t91 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t98 * mrSges(6,2) - t94 * mrSges(6,3) + Ifges(6,1) * t113 / 0.2e1) * t113 + (t117 * mrSges(5,2) - t100 * mrSges(5,3) + t131 * Ifges(5,5) + Ifges(5,1) * t119 / 0.2e1) * t119 + (-t98 * mrSges(6,1) + t95 * mrSges(6,3) + Ifges(6,4) * t113 + Ifges(6,2) * t112 / 0.2e1) * t112 + (t96 * mrSges(7,2) - t90 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t106 / 0.2e1) * t106 + ((-t129 * mrSges(3,3) + Ifges(3,5) * t135 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t144 / 0.2e1) * t151) * t144 + (-t114 * mrSges(4,1) + t115 * mrSges(4,2) + t130 * mrSges(3,3) - Ifges(4,5) * t128 + Ifges(3,6) * t135 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t144 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t146) * t151) * t146) * t151 + (-Ifges(4,6) * t148 - t121 * mrSges(4,1) + t115 * mrSges(4,3) + Ifges(4,4) * t128 + Ifges(4,2) * t127 / 0.2e1) * t127 + (-t96 * mrSges(7,1) + t91 * mrSges(7,3) + Ifges(7,4) * t106 + Ifges(7,6) * t116 + Ifges(7,2) * t105 / 0.2e1) * t105 + (t117 * mrSges(5,1) + t94 * mrSges(6,1) - t95 * mrSges(6,2) - t101 * mrSges(5,3) - Ifges(5,4) * t119 + Ifges(6,5) * t113 - Ifges(5,6) * t131 + Ifges(6,6) * t112 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t118) * t118;
T  = t1;
