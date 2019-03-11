% Calculate kinetic energy for
% S6RRPRRR5
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:37:57
% EndTime: 2019-03-09 13:37:58
% DurationCPUTime: 0.71s
% Computational Cost: add. (1176->115), mult. (3112->183), div. (0->0), fcn. (2514->12), ass. (0->53)
t134 = sin(pkin(12));
t136 = cos(pkin(12));
t141 = sin(qJ(2));
t145 = cos(qJ(2));
t135 = sin(pkin(6));
t152 = qJD(1) * t135;
t125 = (t134 * t141 - t136 * t145) * t152;
t139 = sin(qJ(5));
t143 = cos(qJ(5));
t151 = cos(pkin(6)) * qJD(1);
t150 = pkin(1) * t151;
t132 = t145 * t150;
t133 = qJD(2) + t151;
t149 = t141 * t152;
t119 = pkin(2) * t133 + t132 + (-pkin(8) - qJ(3)) * t149;
t148 = t145 * t152;
t128 = pkin(8) * t148 + t141 * t150;
t122 = qJ(3) * t148 + t128;
t109 = t134 * t119 + t136 * t122;
t107 = pkin(9) * t133 + t109;
t126 = (t134 * t145 + t136 * t141) * t152;
t129 = qJD(3) + (-pkin(2) * t145 - pkin(1)) * t152;
t113 = pkin(3) * t125 - pkin(9) * t126 + t129;
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t101 = t144 * t107 + t140 * t113;
t124 = qJD(4) + t125;
t96 = pkin(10) * t124 + t101;
t108 = t119 * t136 - t134 * t122;
t106 = -pkin(3) * t133 - t108;
t117 = -t140 * t126 + t133 * t144;
t118 = t126 * t144 + t133 * t140;
t99 = -pkin(4) * t117 - pkin(10) * t118 + t106;
t92 = t139 * t99 + t143 * t96;
t91 = -t139 * t96 + t143 * t99;
t100 = -t140 * t107 + t113 * t144;
t116 = qJD(5) - t117;
t95 = -pkin(4) * t124 - t100;
t146 = qJD(1) ^ 2;
t142 = cos(qJ(6));
t138 = sin(qJ(6));
t127 = -pkin(8) * t149 + t132;
t114 = qJD(6) + t116;
t111 = t118 * t143 + t124 * t139;
t110 = -t118 * t139 + t124 * t143;
t103 = t110 * t138 + t111 * t142;
t102 = t110 * t142 - t111 * t138;
t93 = -pkin(5) * t110 + t95;
t90 = pkin(11) * t110 + t92;
t89 = pkin(5) * t116 - pkin(11) * t111 + t91;
t88 = t138 * t89 + t142 * t90;
t87 = -t138 * t90 + t142 * t89;
t1 = m(5) * (t100 ^ 2 + t101 ^ 2 + t106 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t129 ^ 2) / 0.2e1 + t146 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t135 ^ 2 * t146 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + (t129 * mrSges(4,2) - t108 * mrSges(4,3) + Ifges(4,1) * t126 / 0.2e1) * t126 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t124 / 0.2e1) * t124 + (t91 * mrSges(6,1) - t92 * mrSges(6,2) + Ifges(6,3) * t116 / 0.2e1) * t116 + (t87 * mrSges(7,1) - t88 * mrSges(7,2) + Ifges(7,3) * t114 / 0.2e1) * t114 - (-t129 * mrSges(4,1) + t109 * mrSges(4,3) + Ifges(4,4) * t126 - Ifges(4,2) * t125 / 0.2e1) * t125 + (t106 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,5) * t124 + Ifges(5,1) * t118 / 0.2e1) * t118 + (t95 * mrSges(6,2) - t91 * mrSges(6,3) + Ifges(6,5) * t116 + Ifges(6,1) * t111 / 0.2e1) * t111 + (t93 * mrSges(7,2) - t87 * mrSges(7,3) + Ifges(7,5) * t114 + Ifges(7,1) * t103 / 0.2e1) * t103 + (-t106 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t118 + Ifges(5,6) * t124 + Ifges(5,2) * t117 / 0.2e1) * t117 + (-t95 * mrSges(6,1) + t92 * mrSges(6,3) + Ifges(6,4) * t111 + Ifges(6,6) * t116 + Ifges(6,2) * t110 / 0.2e1) * t110 + (-t93 * mrSges(7,1) + t88 * mrSges(7,3) + Ifges(7,4) * t103 + Ifges(7,6) * t114 + Ifges(7,2) * t102 / 0.2e1) * t102 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t145 / 0.2e1) * t145 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t145 + Ifges(3,1) * t141 / 0.2e1) * t141) * t152 + (-t127 * t141 + t128 * t145) * mrSges(3,3)) * t152 + (t127 * mrSges(3,1) + t108 * mrSges(4,1) - t128 * mrSges(3,2) - t109 * mrSges(4,2) + Ifges(4,5) * t126 - Ifges(4,6) * t125 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t133 + (Ifges(3,5) * t141 + Ifges(3,6) * t145) * t152) * t133;
T  = t1;
