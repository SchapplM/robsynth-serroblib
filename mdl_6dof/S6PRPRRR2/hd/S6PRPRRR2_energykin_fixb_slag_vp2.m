% Calculate kinetic energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:35
% EndTime: 2019-03-08 20:26:35
% DurationCPUTime: 0.34s
% Computational Cost: add. (367->84), mult. (789->136), div. (0->0), fcn. (544->12), ass. (0->41)
t131 = cos(qJ(2));
t121 = sin(pkin(6));
t137 = qJD(1) * t121;
t112 = qJD(2) * pkin(2) + t131 * t137;
t120 = sin(pkin(12));
t122 = cos(pkin(12));
t127 = sin(qJ(2));
t134 = t127 * t137;
t107 = t112 * t122 - t120 * t134;
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t100 = (-pkin(4) * t130 - pkin(9) * t126 - pkin(3)) * qJD(2) - t107;
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t108 = t120 * t112 + t122 * t134;
t106 = qJD(2) * pkin(8) + t108;
t123 = cos(pkin(6));
t117 = qJD(1) * t123 + qJD(3);
t97 = t130 * t106 + t126 * t117;
t95 = qJD(4) * pkin(9) + t97;
t91 = t125 * t100 + t129 * t95;
t136 = qJD(2) * t126;
t135 = qJD(2) * t130;
t90 = t129 * t100 - t125 * t95;
t96 = -t126 * t106 + t117 * t130;
t118 = qJD(5) - t135;
t94 = -qJD(4) * pkin(4) - t96;
t128 = cos(qJ(6));
t124 = sin(qJ(6));
t116 = qJD(6) + t118;
t111 = qJD(4) * t125 + t129 * t136;
t110 = qJD(4) * t129 - t125 * t136;
t105 = -qJD(2) * pkin(3) - t107;
t102 = t110 * t124 + t111 * t128;
t101 = t110 * t128 - t111 * t124;
t92 = -pkin(5) * t110 + t94;
t89 = pkin(10) * t110 + t91;
t88 = pkin(5) * t118 - pkin(10) * t111 + t90;
t87 = t124 * t88 + t128 * t89;
t86 = -t124 * t89 + t128 * t88;
t1 = m(4) * (t107 ^ 2 + t108 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t94 ^ 2) / 0.2e1 + (t90 * mrSges(6,1) - t91 * mrSges(6,2) + Ifges(6,3) * t118 / 0.2e1) * t118 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t96 * mrSges(5,1) - t97 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t94 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,5) * t118 + Ifges(6,1) * t111 / 0.2e1) * t111 + (t92 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t102 / 0.2e1) * t102 + (m(3) * (t123 ^ 2 + (t127 ^ 2 + t131 ^ 2) * t121 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t94 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t111 + Ifges(6,6) * t118 + Ifges(6,2) * t110 / 0.2e1) * t110 + (-t92 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,4) * t102 + Ifges(7,6) * t116 + Ifges(7,2) * t101 / 0.2e1) * t101 + (t107 * mrSges(4,1) - t108 * mrSges(4,2) + (mrSges(3,1) * t131 - mrSges(3,2) * t127) * t137 + (-t105 * mrSges(5,1) + t97 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t135 / 0.2e1) * t130 + (t105 * mrSges(5,2) - t96 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t126 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t130 + Ifges(5,1) * t126 / 0.2e1) * t126) * qJD(2)) * qJD(2);
T  = t1;
