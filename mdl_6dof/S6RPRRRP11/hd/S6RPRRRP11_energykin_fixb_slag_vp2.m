% Calculate kinetic energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:22
% EndTime: 2019-03-09 06:34:23
% DurationCPUTime: 0.70s
% Computational Cost: add. (1372->120), mult. (4296->190), div. (0->0), fcn. (3588->12), ass. (0->52)
t137 = sin(pkin(12));
t139 = sin(pkin(6));
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t140 = cos(pkin(12));
t141 = cos(pkin(7));
t156 = t140 * t141;
t138 = sin(pkin(7));
t142 = cos(pkin(6));
t158 = t138 * t142;
t125 = ((-t137 * t145 + t148 * t156) * t139 + t148 * t158) * qJD(1);
t154 = qJD(1) * t139;
t151 = qJ(2) * t154;
t153 = pkin(1) * qJD(1) * t142;
t132 = t137 * t153 + t140 * t151;
t122 = (t139 * t156 + t158) * qJD(1) * pkin(9) + t132;
t135 = t140 * t153;
t124 = t135 + (pkin(2) * t142 + (-pkin(9) * t141 - qJ(2)) * t139 * t137) * qJD(1);
t129 = qJD(2) + (-pkin(9) * t137 * t138 - pkin(2) * t140 - pkin(1)) * t154;
t111 = -t145 * t122 + (t124 * t141 + t129 * t138) * t148;
t130 = qJD(3) + (-t138 * t139 * t140 + t141 * t142) * qJD(1);
t109 = -t130 * pkin(3) - t111;
t155 = t141 * t145;
t157 = t138 * t145;
t126 = (t142 * t157 + (t137 * t148 + t140 * t155) * t139) * qJD(1);
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t117 = -t126 * t144 + t130 * t147;
t118 = t126 * t147 + t130 * t144;
t104 = -t117 * pkin(4) - t118 * pkin(11) + t109;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t115 = -t124 * t138 + t129 * t141;
t106 = -pkin(3) * t125 - pkin(10) * t126 + t115;
t112 = t122 * t148 + t124 * t155 + t129 * t157;
t110 = pkin(10) * t130 + t112;
t101 = t106 * t144 + t110 * t147;
t123 = qJD(4) - t125;
t99 = pkin(11) * t123 + t101;
t95 = t104 * t143 + t146 * t99;
t94 = t104 * t146 - t143 * t99;
t100 = t106 * t147 - t110 * t144;
t98 = -pkin(4) * t123 - t100;
t136 = -pkin(1) * t154 + qJD(2);
t131 = -t137 * t151 + t135;
t116 = qJD(5) - t117;
t114 = t118 * t146 + t123 * t143;
t113 = -t118 * t143 + t123 * t146;
t96 = -pkin(5) * t113 + qJD(6) + t98;
t93 = qJ(6) * t113 + t95;
t92 = pkin(5) * t116 - qJ(6) * t114 + t94;
t1 = m(3) * (t131 ^ 2 + t132 ^ 2 + t136 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t92 ^ 2 + t93 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t109 ^ 2) / 0.2e1 + (t111 * mrSges(4,1) - t112 * mrSges(4,2) + Ifges(4,3) * t130 / 0.2e1) * t130 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t123 / 0.2e1) * t123 + (t115 * mrSges(4,2) - t111 * mrSges(4,3) + Ifges(4,5) * t130 + Ifges(4,1) * t126 / 0.2e1) * t126 + (t109 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,5) * t123 + Ifges(5,1) * t118 / 0.2e1) * t118 + (-t115 * mrSges(4,1) + t112 * mrSges(4,3) + Ifges(4,4) * t126 + Ifges(4,6) * t130 + Ifges(4,2) * t125 / 0.2e1) * t125 + (-t109 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t118 + Ifges(5,6) * t123 + Ifges(5,2) * t117 / 0.2e1) * t117 + (t94 * mrSges(6,1) + t92 * mrSges(7,1) - t95 * mrSges(6,2) - t93 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t116) * t116 + (t98 * mrSges(6,2) + t96 * mrSges(7,2) - t94 * mrSges(6,3) - t92 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t114 + (Ifges(6,5) + Ifges(7,5)) * t116) * t114 + (-t98 * mrSges(6,1) - t96 * mrSges(7,1) + t95 * mrSges(6,3) + t93 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t113 + (Ifges(6,6) + Ifges(7,6)) * t116 + (Ifges(6,4) + Ifges(7,4)) * t114) * t113 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t136 * (-mrSges(3,1) * t140 + mrSges(3,2) * t137) + (Ifges(3,2) * t140 ^ 2 / 0.2e1 + (Ifges(3,4) * t140 + Ifges(3,1) * t137 / 0.2e1) * t137) * t154 + (-t131 * t137 + t132 * t140) * mrSges(3,3)) * t139 + (t131 * mrSges(3,1) - t132 * mrSges(3,2) + (Ifges(3,3) * t142 / 0.2e1 + (Ifges(3,5) * t137 + Ifges(3,6) * t140) * t139) * qJD(1)) * t142) * qJD(1);
T  = t1;
