% Calculate kinetic energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:07
% EndTime: 2019-03-09 01:39:07
% DurationCPUTime: 0.46s
% Computational Cost: add. (501->93), mult. (1154->147), div. (0->0), fcn. (796->10), ass. (0->41)
t125 = m(3) / 0.2e1;
t124 = cos(qJ(4));
t111 = sin(pkin(11));
t114 = cos(pkin(11));
t118 = sin(qJ(4));
t113 = sin(pkin(9));
t107 = (pkin(1) * t113 + qJ(3)) * qJD(1);
t115 = cos(pkin(10));
t110 = t115 * qJD(2);
t112 = sin(pkin(10));
t97 = t110 + (-pkin(7) * qJD(1) - t107) * t112;
t100 = t112 * qJD(2) + t115 * t107;
t123 = t115 * qJD(1);
t98 = pkin(7) * t123 + t100;
t88 = t118 * t97 + t124 * t98;
t84 = qJD(4) * qJ(5) + t88;
t116 = cos(pkin(9));
t122 = -pkin(1) * t116 - pkin(2);
t102 = qJD(3) + (-pkin(3) * t115 + t122) * qJD(1);
t103 = qJD(1) * t112 * t118 - t124 * t123;
t104 = (t112 * t124 + t115 * t118) * qJD(1);
t91 = pkin(4) * t103 - qJ(5) * t104 + t102;
t80 = t111 * t91 + t114 * t84;
t79 = -t111 * t84 + t114 * t91;
t87 = -t118 * t98 + t124 * t97;
t83 = -qJD(4) * pkin(4) + qJD(5) - t87;
t119 = cos(qJ(6));
t117 = sin(qJ(6));
t106 = qJD(1) * t122 + qJD(3);
t101 = qJD(6) + t103;
t99 = -t107 * t112 + t110;
t96 = qJD(4) * t111 + t104 * t114;
t95 = qJD(4) * t114 - t104 * t111;
t86 = t117 * t95 + t119 * t96;
t85 = -t117 * t96 + t119 * t95;
t81 = -t95 * pkin(5) + t83;
t78 = pkin(8) * t95 + t80;
t77 = pkin(5) * t103 - pkin(8) * t96 + t79;
t76 = t117 * t77 + t119 * t78;
t75 = -t117 * t78 + t119 * t77;
t1 = qJD(2) ^ 2 * t125 + m(5) * (t102 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t106 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + (t83 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t96 / 0.2e1) * t96 + (t81 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,1) * t86 / 0.2e1) * t86 + (t102 * mrSges(5,2) - t87 * mrSges(5,3) + Ifges(5,1) * t104 / 0.2e1) * t104 + (-t83 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t96 + Ifges(6,2) * t95 / 0.2e1) * t95 + (-t81 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t87 * mrSges(5,1) - t88 * mrSges(5,2) + Ifges(5,5) * t104 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + t86 * Ifges(7,5) + t85 * Ifges(7,6) + Ifges(7,3) * t101 / 0.2e1) * t101 + (t102 * mrSges(5,1) + t79 * mrSges(6,1) - t80 * mrSges(6,2) - t88 * mrSges(5,3) - Ifges(5,4) * t104 + Ifges(6,5) * t96 - Ifges(5,6) * qJD(4) + Ifges(6,6) * t95 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t103) * t103 + (t106 * (-mrSges(4,1) * t115 + mrSges(4,2) * t112) + (t100 * t115 - t99 * t112) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t116 * mrSges(3,1) - t113 * mrSges(3,2) + (t113 ^ 2 + t116 ^ 2) * t125 * pkin(1)) * pkin(1) + Ifges(4,2) * t115 ^ 2 / 0.2e1 + (Ifges(4,4) * t115 + Ifges(4,1) * t112 / 0.2e1) * t112) * qJD(1)) * qJD(1);
T  = t1;
