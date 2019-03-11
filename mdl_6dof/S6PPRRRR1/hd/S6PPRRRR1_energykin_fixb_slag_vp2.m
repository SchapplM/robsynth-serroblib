% Calculate kinetic energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:39
% EndTime: 2019-03-08 18:59:39
% DurationCPUTime: 0.36s
% Computational Cost: add. (422->83), mult. (973->138), div. (0->0), fcn. (780->14), ass. (0->45)
t128 = cos(pkin(6)) * qJD(1) + qJD(2);
t132 = sin(pkin(7));
t135 = cos(pkin(7));
t134 = cos(pkin(13));
t133 = sin(pkin(6));
t152 = qJD(1) * t133;
t149 = t134 * t152;
t154 = t128 * t132 + t135 * t149;
t137 = sin(qJ(5));
t138 = sin(qJ(4));
t141 = cos(qJ(5));
t142 = cos(qJ(4));
t122 = (t137 * t138 - t141 * t142) * qJD(3);
t139 = sin(qJ(3));
t143 = cos(qJ(3));
t131 = sin(pkin(13));
t150 = t131 * t152;
t114 = -t139 * t150 + t154 * t143;
t115 = t154 * t139 + t143 * t150;
t113 = qJD(3) * pkin(9) + t115;
t120 = t135 * t128 - t132 * t149;
t119 = t142 * t120;
t106 = qJD(4) * pkin(4) + t119 + (-pkin(10) * qJD(3) - t113) * t138;
t109 = t142 * t113 + t138 * t120;
t151 = qJD(3) * t142;
t107 = pkin(10) * t151 + t109;
t102 = t137 * t106 + t141 * t107;
t101 = t141 * t106 - t137 * t107;
t110 = (-pkin(4) * t142 - pkin(3)) * qJD(3) - t114;
t144 = qJD(1) ^ 2;
t140 = cos(qJ(6));
t136 = sin(qJ(6));
t130 = qJD(4) + qJD(5);
t123 = (t137 * t142 + t138 * t141) * qJD(3);
t121 = qJD(6) + t122;
t117 = t140 * t123 + t136 * t130;
t116 = -t136 * t123 + t140 * t130;
t112 = -qJD(3) * pkin(3) - t114;
t108 = -t138 * t113 + t119;
t104 = t122 * pkin(5) - t123 * pkin(11) + t110;
t100 = t130 * pkin(11) + t102;
t99 = -t130 * pkin(5) - t101;
t98 = t140 * t100 + t136 * t104;
t97 = -t136 * t100 + t140 * t104;
t1 = m(3) * (t128 ^ 2 + (t131 ^ 2 + t134 ^ 2) * t144 * t133 ^ 2) / 0.2e1 + m(2) * t144 / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t120 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + (t101 * mrSges(6,1) - t102 * mrSges(6,2) + Ifges(6,3) * t130 / 0.2e1) * t130 + (t97 * mrSges(7,1) - t98 * mrSges(7,2) + Ifges(7,3) * t121 / 0.2e1) * t121 + (t108 * mrSges(5,1) - t109 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t110 * mrSges(6,2) - t101 * mrSges(6,3) + Ifges(6,5) * t130 + Ifges(6,1) * t123 / 0.2e1) * t123 + (t99 * mrSges(7,2) - t97 * mrSges(7,3) + Ifges(7,5) * t121 + Ifges(7,1) * t117 / 0.2e1) * t117 - (-t110 * mrSges(6,1) + t102 * mrSges(6,3) + Ifges(6,4) * t123 + Ifges(6,6) * t130 - Ifges(6,2) * t122 / 0.2e1) * t122 + (-t99 * mrSges(7,1) + t98 * mrSges(7,3) + Ifges(7,4) * t117 + Ifges(7,6) * t121 + Ifges(7,2) * t116 / 0.2e1) * t116 + (t114 * mrSges(4,1) - t115 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t112 * mrSges(5,1) + t109 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t151 / 0.2e1) * t142 + (t112 * mrSges(5,2) - t108 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t142 + Ifges(5,1) * t138 / 0.2e1) * qJD(3)) * t138) * qJD(3);
T  = t1;
