% Calculate kinetic energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:59
% EndTime: 2019-12-31 21:37:00
% DurationCPUTime: 0.40s
% Computational Cost: add. (646->90), mult. (1479->146), div. (0->0), fcn. (1132->10), ass. (0->42)
t126 = cos(qJ(3));
t111 = sin(pkin(10));
t113 = cos(pkin(10));
t124 = cos(pkin(5)) * qJD(1);
t110 = qJD(2) + t124;
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t112 = sin(pkin(5));
t125 = qJD(1) * t112;
t122 = t117 * t125;
t101 = -t126 * t110 + t116 * t122;
t102 = t116 * t110 + t126 * t122;
t119 = cos(qJ(2));
t123 = pkin(1) * t124;
t103 = -pkin(7) * t122 + t119 * t123;
t97 = -pkin(2) * t110 - t103;
t86 = pkin(3) * t101 - qJ(4) * t102 + t97;
t121 = t119 * t125;
t106 = qJD(3) - t121;
t104 = pkin(7) * t121 + t117 * t123;
t98 = pkin(8) * t110 + t104;
t99 = (-pkin(2) * t119 - pkin(8) * t117 - pkin(1)) * t125;
t91 = t116 * t99 + t126 * t98;
t89 = qJ(4) * t106 + t91;
t80 = t111 * t86 + t113 * t89;
t79 = -t111 * t89 + t113 * t86;
t90 = -t116 * t98 + t126 * t99;
t88 = -t106 * pkin(3) + qJD(4) - t90;
t120 = qJD(1) ^ 2;
t118 = cos(qJ(5));
t115 = sin(qJ(5));
t100 = qJD(5) + t101;
t93 = t102 * t113 + t106 * t111;
t92 = -t102 * t111 + t106 * t113;
t83 = t115 * t92 + t118 * t93;
t82 = -t115 * t93 + t118 * t92;
t81 = -t92 * pkin(4) + t88;
t78 = pkin(9) * t92 + t80;
t77 = pkin(4) * t101 - pkin(9) * t93 + t79;
t76 = t115 * t77 + t118 * t78;
t75 = -t115 * t78 + t118 * t77;
t1 = m(3) * (pkin(1) ^ 2 * t112 ^ 2 * t120 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t120 * Ifges(2,3) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + (t88 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,1) * t83 / 0.2e1) * t83 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,3) * t106 / 0.2e1) * t106 + (-t88 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t83 + Ifges(6,2) * t82 / 0.2e1) * t82 + (t97 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,5) * t106 + Ifges(4,1) * t102 / 0.2e1) * t102 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t119 / 0.2e1) * t119 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t119 + Ifges(3,1) * t117 / 0.2e1) * t117) * t125 + (-t103 * t117 + t104 * t119) * mrSges(3,3)) * t125 + (t103 * mrSges(3,1) - t104 * mrSges(3,2) + Ifges(3,3) * t110 / 0.2e1 + (Ifges(3,5) * t117 + Ifges(3,6) * t119) * t125) * t110 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,5) * t83 + Ifges(6,6) * t82 + Ifges(6,3) * t100 / 0.2e1) * t100 + (t97 * mrSges(4,1) + t79 * mrSges(5,1) - t80 * mrSges(5,2) - t91 * mrSges(4,3) - Ifges(4,4) * t102 + Ifges(5,5) * t93 - Ifges(4,6) * t106 + Ifges(5,6) * t92 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101) * t101;
T = t1;
