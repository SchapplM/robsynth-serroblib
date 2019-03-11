% Calculate kinetic energy for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:02
% EndTime: 2019-03-09 00:13:02
% DurationCPUTime: 0.41s
% Computational Cost: add. (454->94), mult. (951->141), div. (0->0), fcn. (661->10), ass. (0->39)
t129 = cos(qJ(5));
t116 = sin(qJ(5));
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t118 = sin(qJ(3));
t126 = qJD(2) * t118;
t107 = qJD(3) * t117 + t120 * t126;
t121 = cos(qJ(3));
t125 = qJD(2) * t121;
t112 = qJD(4) - t125;
t122 = cos(qJ(2));
t114 = sin(pkin(6));
t128 = qJD(1) * t114;
t124 = t122 * t128;
t102 = -t124 + (-pkin(3) * t121 - pkin(9) * t118 - pkin(2)) * qJD(2);
t119 = sin(qJ(2));
t108 = qJD(2) * pkin(8) + t119 * t128;
t115 = cos(pkin(6));
t127 = qJD(1) * t115;
t101 = t121 * t108 + t118 * t127;
t97 = qJD(3) * pkin(9) + t101;
t90 = t120 * t102 - t117 * t97;
t87 = pkin(4) * t112 - pkin(10) * t107 + t90;
t106 = qJD(3) * t120 - t117 * t126;
t91 = t117 * t102 + t120 * t97;
t89 = pkin(10) * t106 + t91;
t84 = t116 * t87 + t129 * t89;
t100 = -t118 * t108 + t121 * t127;
t83 = -t116 * t89 + t129 * t87;
t96 = -qJD(3) * pkin(3) - t100;
t92 = -pkin(4) * t106 + t96;
t110 = qJD(5) + t112;
t109 = -qJD(2) * pkin(2) - t124;
t94 = t116 * t106 + t129 * t107;
t93 = -t129 * t106 + t107 * t116;
t85 = pkin(5) * t93 - qJ(6) * t94 + t92;
t82 = qJ(6) * t110 + t84;
t81 = -t110 * pkin(5) + qJD(6) - t83;
t1 = m(4) * (t100 ^ 2 + t101 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t112 / 0.2e1) * t112 + (t100 * mrSges(4,1) - t101 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t96 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t112 + Ifges(5,1) * t107 / 0.2e1) * t107 + (m(3) * (t115 ^ 2 + (t119 ^ 2 + t122 ^ 2) * t114 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t96 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t107 + Ifges(5,6) * t112 + Ifges(5,2) * t106 / 0.2e1) * t106 + (t92 * mrSges(6,2) + t81 * mrSges(7,2) - t83 * mrSges(6,3) - t85 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t94) * t94 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t122 - mrSges(3,2) * t119) * t128 + (-t109 * mrSges(4,1) + t101 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t125 / 0.2e1) * t121 + (t109 * mrSges(4,2) - t100 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t121 + Ifges(4,1) * t118 / 0.2e1) * qJD(2)) * t118) * qJD(2) + (t92 * mrSges(6,1) + t85 * mrSges(7,1) - t82 * mrSges(7,2) - t84 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t93 + (-Ifges(6,4) + Ifges(7,5)) * t94) * t93 + (t83 * mrSges(6,1) - t81 * mrSges(7,1) - t84 * mrSges(6,2) + t82 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t110 + (Ifges(7,4) + Ifges(6,5)) * t94 + (-Ifges(6,6) + Ifges(7,6)) * t93) * t110;
T  = t1;
