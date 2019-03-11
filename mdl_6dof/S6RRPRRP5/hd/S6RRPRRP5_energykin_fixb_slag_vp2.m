% Calculate kinetic energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:42
% EndTime: 2019-03-09 11:57:42
% DurationCPUTime: 0.63s
% Computational Cost: add. (898->111), mult. (2392->166), div. (0->0), fcn. (1892->10), ass. (0->46)
t120 = sin(pkin(11));
t122 = cos(pkin(11));
t126 = sin(qJ(2));
t129 = cos(qJ(2));
t121 = sin(pkin(6));
t136 = qJD(1) * t121;
t111 = (t120 * t126 - t122 * t129) * t136;
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t110 = qJD(4) + t111;
t112 = (t120 * t129 + t122 * t126) * t136;
t115 = qJD(3) + (-pkin(2) * t129 - pkin(1)) * t136;
t101 = pkin(3) * t111 - pkin(9) * t112 + t115;
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t135 = cos(pkin(6)) * qJD(1);
t119 = qJD(2) + t135;
t134 = pkin(1) * t135;
t118 = t129 * t134;
t133 = t126 * t136;
t106 = pkin(2) * t119 + t118 + (-pkin(8) - qJ(3)) * t133;
t132 = t129 * t136;
t114 = pkin(8) * t132 + t126 * t134;
t109 = qJ(3) * t132 + t114;
t97 = t120 * t106 + t122 * t109;
t95 = pkin(9) * t119 + t97;
t91 = t125 * t101 + t128 * t95;
t86 = pkin(10) * t110 + t91;
t104 = -t112 * t125 + t119 * t128;
t105 = t112 * t128 + t119 * t125;
t96 = t106 * t122 - t120 * t109;
t94 = -pkin(3) * t119 - t96;
t89 = -pkin(4) * t104 - pkin(10) * t105 + t94;
t82 = t124 * t89 + t127 * t86;
t81 = -t124 * t86 + t127 * t89;
t90 = t101 * t128 - t125 * t95;
t85 = -pkin(4) * t110 - t90;
t130 = qJD(1) ^ 2;
t113 = -pkin(8) * t133 + t118;
t103 = qJD(5) - t104;
t99 = t105 * t127 + t110 * t124;
t98 = -t105 * t124 + t110 * t127;
t83 = -pkin(5) * t98 + qJD(6) + t85;
t80 = qJ(6) * t98 + t82;
t79 = pkin(5) * t103 - qJ(6) * t99 + t81;
t1 = m(5) * (t90 ^ 2 + t91 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t121 ^ 2 * t130 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t130 * Ifges(2,3) / 0.2e1 + m(4) * (t115 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + (t115 * mrSges(4,2) - t96 * mrSges(4,3) + Ifges(4,1) * t112 / 0.2e1) * t112 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t110 / 0.2e1) * t110 - (-t115 * mrSges(4,1) + t97 * mrSges(4,3) + Ifges(4,4) * t112 - Ifges(4,2) * t111 / 0.2e1) * t111 + (t94 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t110 + Ifges(5,1) * t105 / 0.2e1) * t105 + (-t94 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,6) * t110 + Ifges(5,2) * t104 / 0.2e1) * t104 + (t85 * mrSges(6,2) + t83 * mrSges(7,2) - t81 * mrSges(6,3) - t79 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t99) * t99 + (-t85 * mrSges(6,1) - t83 * mrSges(7,1) + t82 * mrSges(6,3) + t80 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t98 + (Ifges(6,4) + Ifges(7,4)) * t99) * t98 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t129 / 0.2e1) * t129 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t129 + Ifges(3,1) * t126 / 0.2e1) * t126) * t136 + (-t113 * t126 + t114 * t129) * mrSges(3,3)) * t136 + (t113 * mrSges(3,1) + t96 * mrSges(4,1) - t114 * mrSges(3,2) - t97 * mrSges(4,2) + Ifges(4,5) * t112 - Ifges(4,6) * t111 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t119 + (Ifges(3,5) * t126 + Ifges(3,6) * t129) * t136) * t119 + (t81 * mrSges(6,1) + t79 * mrSges(7,1) - t82 * mrSges(6,2) - t80 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t103 + (Ifges(6,5) + Ifges(7,5)) * t99 + (Ifges(6,6) + Ifges(7,6)) * t98) * t103;
T  = t1;
