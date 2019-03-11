% Calculate kinetic energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:33
% EndTime: 2019-03-10 02:59:33
% DurationCPUTime: 0.69s
% Computational Cost: add. (1666->119), mult. (4274->184), div. (0->0), fcn. (3564->12), ass. (0->54)
t144 = sin(qJ(2));
t147 = cos(qJ(2));
t138 = sin(pkin(6));
t156 = qJD(1) * t138;
t151 = t147 * t156;
t155 = cos(pkin(6)) * qJD(1);
t154 = pkin(1) * t155;
t132 = pkin(9) * t151 + t144 * t154;
t139 = cos(pkin(7));
t136 = qJD(2) + t155;
t137 = sin(pkin(7));
t159 = t136 * t137;
t122 = (t139 * t151 + t159) * pkin(10) + t132;
t129 = (-pkin(10) * t137 * t144 - pkin(2) * t147 - pkin(1)) * t156;
t143 = sin(qJ(3));
t146 = cos(qJ(3));
t135 = t147 * t154;
t152 = t144 * t156;
t124 = pkin(2) * t136 + t135 + (-pkin(10) * t139 - pkin(9)) * t152;
t160 = t124 * t139;
t110 = -t143 * t122 + (t129 * t137 + t160) * t146;
t157 = t139 * t147;
t125 = (-t143 * t144 + t146 * t157) * t156 + t146 * t159;
t161 = cos(qJ(5));
t130 = t136 * t139 - t137 * t151 + qJD(3);
t108 = -pkin(3) * t130 - t110;
t158 = t137 * t143;
t126 = t136 * t158 + (t143 * t157 + t144 * t146) * t156;
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t116 = -t126 * t142 + t130 * t145;
t117 = t126 * t145 + t130 * t142;
t103 = -pkin(4) * t116 - pkin(12) * t117 + t108;
t141 = sin(qJ(5));
t114 = -t124 * t137 + t129 * t139;
t105 = -pkin(3) * t125 - pkin(11) * t126 + t114;
t111 = t122 * t146 + t129 * t158 + t143 * t160;
t109 = pkin(11) * t130 + t111;
t101 = t105 * t142 + t109 * t145;
t123 = qJD(4) - t125;
t99 = pkin(12) * t123 + t101;
t95 = t103 * t141 + t161 * t99;
t100 = t105 * t145 - t109 * t142;
t98 = -pkin(4) * t123 - t100;
t94 = t103 * t161 - t141 * t99;
t148 = qJD(1) ^ 2;
t131 = -pkin(9) * t152 + t135;
t115 = qJD(5) - t116;
t113 = t117 * t161 + t123 * t141;
t112 = t117 * t141 - t123 * t161;
t96 = pkin(5) * t112 - qJ(6) * t113 + t98;
t93 = qJ(6) * t115 + t95;
t92 = -pkin(5) * t115 + qJD(6) - t94;
t1 = t148 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t138 ^ 2 * t148 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t92 ^ 2 + t93 ^ 2 + t96 ^ 2) / 0.2e1 + (t110 * mrSges(4,1) - t111 * mrSges(4,2) + Ifges(4,3) * t130 / 0.2e1) * t130 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t123 / 0.2e1) * t123 + (t114 * mrSges(4,2) - t110 * mrSges(4,3) + Ifges(4,5) * t130 + Ifges(4,1) * t126 / 0.2e1) * t126 + (t108 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,5) * t123 + Ifges(5,1) * t117 / 0.2e1) * t117 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t147 / 0.2e1) * t147 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t147 + Ifges(3,1) * t144 / 0.2e1) * t144) * t156 + (-t131 * t144 + t132 * t147) * mrSges(3,3)) * t156 + (t131 * mrSges(3,1) - t132 * mrSges(3,2) + Ifges(3,3) * t136 / 0.2e1 + (Ifges(3,5) * t144 + Ifges(3,6) * t147) * t156) * t136 + (-t114 * mrSges(4,1) + t111 * mrSges(4,3) + Ifges(4,4) * t126 + Ifges(4,6) * t130 + Ifges(4,2) * t125 / 0.2e1) * t125 + (-t108 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t117 + Ifges(5,6) * t123 + Ifges(5,2) * t116 / 0.2e1) * t116 + (t94 * mrSges(6,1) - t92 * mrSges(7,1) - t95 * mrSges(6,2) + t93 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t115) * t115 + (t98 * mrSges(6,2) + t92 * mrSges(7,2) - t94 * mrSges(6,3) - t96 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t113 + (Ifges(7,4) + Ifges(6,5)) * t115) * t113 + (t98 * mrSges(6,1) + t96 * mrSges(7,1) - t93 * mrSges(7,2) - t95 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t112 + (-Ifges(6,6) + Ifges(7,6)) * t115 + (-Ifges(6,4) + Ifges(7,5)) * t113) * t112;
T  = t1;
