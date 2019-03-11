% Calculate kinetic energy for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:10
% EndTime: 2019-03-09 07:12:11
% DurationCPUTime: 0.74s
% Computational Cost: add. (927->104), mult. (2181->164), div. (0->0), fcn. (1674->10), ass. (0->48)
t122 = cos(pkin(11));
t138 = t122 ^ 2;
t137 = m(3) / 0.2e1;
t136 = pkin(7) + qJ(2);
t124 = sin(qJ(5));
t128 = cos(qJ(5));
t121 = sin(pkin(11));
t126 = sin(qJ(3));
t130 = cos(qJ(3));
t113 = (t121 * t130 + t122 * t126) * qJD(1);
t125 = sin(qJ(4));
t129 = cos(qJ(4));
t105 = qJD(3) * t125 + t113 * t129;
t133 = qJD(1) * t122;
t134 = qJD(1) * t121;
t112 = -t126 * t134 + t130 * t133;
t108 = qJD(4) - t112;
t114 = t136 * t134;
t115 = t136 * t133;
t103 = -t126 * t114 + t130 * t115;
t101 = qJD(3) * pkin(8) + t103;
t116 = qJD(2) + (-pkin(2) * t122 - pkin(1)) * qJD(1);
t98 = -pkin(3) * t112 - pkin(8) * t113 + t116;
t91 = -t101 * t125 + t129 * t98;
t85 = pkin(4) * t108 - pkin(9) * t105 + t91;
t104 = qJD(3) * t129 - t113 * t125;
t92 = t129 * t101 + t125 * t98;
t90 = pkin(9) * t104 + t92;
t82 = t124 * t85 + t128 * t90;
t81 = -t124 * t90 + t128 * t85;
t102 = -t114 * t130 - t126 * t115;
t100 = -qJD(3) * pkin(3) - t102;
t107 = qJD(5) + t108;
t93 = -pkin(4) * t104 + t100;
t127 = cos(qJ(6));
t123 = sin(qJ(6));
t118 = -qJD(1) * pkin(1) + qJD(2);
t106 = qJD(6) + t107;
t95 = t104 * t124 + t105 * t128;
t94 = t104 * t128 - t105 * t124;
t88 = -pkin(5) * t94 + t93;
t87 = t123 * t94 + t127 * t95;
t86 = -t123 * t95 + t127 * t94;
t80 = pkin(10) * t94 + t82;
t79 = pkin(5) * t107 - pkin(10) * t95 + t81;
t78 = t123 * t79 + t127 * t80;
t77 = -t123 * t80 + t127 * t79;
t1 = t118 ^ 2 * t137 + m(4) * (t102 ^ 2 + t103 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t88 ^ 2) / 0.2e1 + (t93 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,1) * t95 / 0.2e1) * t95 + (t88 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t87 / 0.2e1) * t87 + (t116 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,1) * t113 / 0.2e1) * t113 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t108 / 0.2e1) * t108 + (-t93 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t95 + Ifges(6,2) * t94 / 0.2e1) * t94 + (-t88 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t87 + Ifges(7,2) * t86 / 0.2e1) * t86 + (-t116 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t113 + Ifges(4,2) * t112 / 0.2e1) * t112 + (t100 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t108 + Ifges(5,1) * t105 / 0.2e1) * t105 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,5) * t95 + Ifges(6,6) * t94 + Ifges(6,3) * t107 / 0.2e1) * t107 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t106 / 0.2e1) * t106 + (-t100 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,6) * t108 + Ifges(5,2) * t104 / 0.2e1) * t104 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,5) * t113 + Ifges(4,6) * t112 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t118 * (-mrSges(3,1) * t122 + mrSges(3,2) * t121) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t137 + mrSges(3,3)) * (t121 ^ 2 + t138) * qJ(2) + Ifges(3,2) * t138 / 0.2e1 + (Ifges(3,4) * t122 + Ifges(3,1) * t121 / 0.2e1) * t121) * qJD(1)) * qJD(1);
T  = t1;
