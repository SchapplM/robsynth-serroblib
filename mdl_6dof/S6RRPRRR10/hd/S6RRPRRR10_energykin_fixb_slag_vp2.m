% Calculate kinetic energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:05
% EndTime: 2019-03-09 14:18:05
% DurationCPUTime: 0.69s
% Computational Cost: add. (1362->114), mult. (3138->180), div. (0->0), fcn. (2538->12), ass. (0->52)
t145 = sin(qJ(2));
t149 = cos(qJ(2));
t139 = sin(pkin(6));
t155 = qJD(1) * t139;
t151 = t149 * t155;
t154 = cos(pkin(6)) * qJD(1);
t153 = pkin(1) * t154;
t132 = pkin(8) * t151 + t145 * t153;
t137 = qJD(2) + t154;
t126 = qJ(3) * t137 + t132;
t128 = (-pkin(2) * t149 - qJ(3) * t145 - pkin(1)) * t155;
t138 = sin(pkin(12));
t140 = cos(pkin(12));
t115 = -t126 * t138 + t140 * t128;
t152 = t145 * t155;
t130 = t137 * t138 + t140 * t152;
t109 = -pkin(3) * t151 - pkin(9) * t130 + t115;
t116 = t140 * t126 + t138 * t128;
t129 = t137 * t140 - t138 * t152;
t112 = pkin(9) * t129 + t116;
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t102 = t144 * t109 + t148 * t112;
t133 = qJD(4) - t151;
t100 = pkin(10) * t133 + t102;
t131 = -pkin(8) * t152 + t149 * t153;
t123 = -pkin(2) * t137 + qJD(3) - t131;
t119 = -pkin(3) * t129 + t123;
t120 = t129 * t148 - t144 * t130;
t121 = t129 * t144 + t130 * t148;
t105 = -pkin(4) * t120 - pkin(10) * t121 + t119;
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t96 = t147 * t100 + t143 * t105;
t95 = -t100 * t143 + t147 * t105;
t101 = t109 * t148 - t144 * t112;
t118 = qJD(5) - t120;
t99 = -pkin(4) * t133 - t101;
t150 = qJD(1) ^ 2;
t146 = cos(qJ(6));
t142 = sin(qJ(6));
t117 = qJD(6) + t118;
t114 = t121 * t147 + t133 * t143;
t113 = -t121 * t143 + t133 * t147;
t107 = t113 * t142 + t114 * t146;
t106 = t113 * t146 - t114 * t142;
t97 = -pkin(5) * t113 + t99;
t94 = pkin(11) * t113 + t96;
t93 = pkin(5) * t118 - pkin(11) * t114 + t95;
t92 = t142 * t93 + t146 * t94;
t91 = -t142 * t94 + t146 * t93;
t1 = m(7) * (t91 ^ 2 + t92 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t95 ^ 2 + t96 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t102 ^ 2 + t119 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t123 ^ 2) / 0.2e1 + t150 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t139 ^ 2 * t150 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + (t131 * mrSges(3,1) - t132 * mrSges(3,2) + Ifges(3,3) * t137 / 0.2e1) * t137 + (t101 * mrSges(5,1) - t102 * mrSges(5,2) + Ifges(5,3) * t133 / 0.2e1) * t133 + (t123 * mrSges(4,2) - t115 * mrSges(4,3) + Ifges(4,1) * t130 / 0.2e1) * t130 + (t95 * mrSges(6,1) - t96 * mrSges(6,2) + Ifges(6,3) * t118 / 0.2e1) * t118 + (t91 * mrSges(7,1) - t92 * mrSges(7,2) + Ifges(7,3) * t117 / 0.2e1) * t117 + (t119 * mrSges(5,2) - t101 * mrSges(5,3) + Ifges(5,5) * t133 + Ifges(5,1) * t121 / 0.2e1) * t121 + (t99 * mrSges(6,2) - t95 * mrSges(6,3) + Ifges(6,5) * t118 + Ifges(6,1) * t114 / 0.2e1) * t114 + (t97 * mrSges(7,2) - t91 * mrSges(7,3) + Ifges(7,5) * t117 + Ifges(7,1) * t107 / 0.2e1) * t107 + ((-t131 * mrSges(3,3) + Ifges(3,5) * t137 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t145 / 0.2e1) * t155) * t145 + (-t115 * mrSges(4,1) + t116 * mrSges(4,2) + t132 * mrSges(3,3) - Ifges(4,5) * t130 + Ifges(3,6) * t137 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t145 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t149) * t155) * t149) * t155 + (-Ifges(4,6) * t151 - t123 * mrSges(4,1) + t116 * mrSges(4,3) + Ifges(4,4) * t130 + Ifges(4,2) * t129 / 0.2e1) * t129 + (-t119 * mrSges(5,1) + t102 * mrSges(5,3) + Ifges(5,4) * t121 + Ifges(5,6) * t133 + Ifges(5,2) * t120 / 0.2e1) * t120 + (-t99 * mrSges(6,1) + t96 * mrSges(6,3) + Ifges(6,4) * t114 + Ifges(6,6) * t118 + Ifges(6,2) * t113 / 0.2e1) * t113 + (-t97 * mrSges(7,1) + t92 * mrSges(7,3) + Ifges(7,4) * t107 + Ifges(7,6) * t117 + Ifges(7,2) * t106 / 0.2e1) * t106;
T  = t1;
