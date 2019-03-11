% Calculate kinetic energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:16
% EndTime: 2019-03-09 00:06:16
% DurationCPUTime: 0.40s
% Computational Cost: add. (456->94), mult. (965->141), div. (0->0), fcn. (675->10), ass. (0->39)
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t118 = sin(qJ(3));
t127 = qJD(2) * t118;
t107 = qJD(3) * t117 + t121 * t127;
t122 = cos(qJ(3));
t126 = qJD(2) * t122;
t112 = qJD(4) - t126;
t123 = cos(qJ(2));
t114 = sin(pkin(6));
t129 = qJD(1) * t114;
t125 = t123 * t129;
t103 = -t125 + (-pkin(3) * t122 - pkin(9) * t118 - pkin(2)) * qJD(2);
t119 = sin(qJ(2));
t108 = qJD(2) * pkin(8) + t119 * t129;
t115 = cos(pkin(6));
t128 = qJD(1) * t115;
t102 = t122 * t108 + t118 * t128;
t98 = qJD(3) * pkin(9) + t102;
t91 = t121 * t103 - t117 * t98;
t87 = pkin(4) * t112 - pkin(10) * t107 + t91;
t106 = qJD(3) * t121 - t117 * t127;
t92 = t117 * t103 + t121 * t98;
t90 = pkin(10) * t106 + t92;
t84 = t116 * t87 + t120 * t90;
t83 = -t116 * t90 + t120 * t87;
t101 = -t118 * t108 + t122 * t128;
t97 = -qJD(3) * pkin(3) - t101;
t93 = -pkin(4) * t106 + t97;
t110 = qJD(5) + t112;
t109 = -qJD(2) * pkin(2) - t125;
t95 = t106 * t116 + t107 * t120;
t94 = t106 * t120 - t107 * t116;
t88 = -pkin(5) * t94 + qJD(6) + t93;
t82 = qJ(6) * t94 + t84;
t81 = pkin(5) * t110 - qJ(6) * t95 + t83;
t1 = m(7) * (t81 ^ 2 + t82 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t93 ^ 2) / 0.2e1 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t112 / 0.2e1) * t112 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t97 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t112 + Ifges(5,1) * t107 / 0.2e1) * t107 + (m(2) / 0.2e1 + m(3) * (t115 ^ 2 + (t119 ^ 2 + t123 ^ 2) * t114 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t97 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t107 + Ifges(5,6) * t112 + Ifges(5,2) * t106 / 0.2e1) * t106 + (t93 * mrSges(6,2) + t88 * mrSges(7,2) - t83 * mrSges(6,3) - t81 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t95) * t95 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t123 - mrSges(3,2) * t119) * t129 + (-t109 * mrSges(4,1) + t102 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t126 / 0.2e1) * t122 + (t109 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t122 + Ifges(4,1) * t118 / 0.2e1) * qJD(2)) * t118) * qJD(2) + (-t93 * mrSges(6,1) - t88 * mrSges(7,1) + t84 * mrSges(6,3) + t82 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t94 + (Ifges(6,4) + Ifges(7,4)) * t95) * t94 + (t83 * mrSges(6,1) + t81 * mrSges(7,1) - t84 * mrSges(6,2) - t82 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t110 + (Ifges(6,5) + Ifges(7,5)) * t95 + (Ifges(6,6) + Ifges(7,6)) * t94) * t110;
T  = t1;
