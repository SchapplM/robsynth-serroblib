% Calculate kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR15_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:53
% EndTime: 2024-09-27 22:23:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (840->96), mult. (2397->160), div. (0->0), fcn. (1940->12), ass. (0->50)
t125 = sin(qJ(3));
t126 = sin(qJ(2));
t129 = cos(qJ(3));
t130 = cos(qJ(2));
t120 = sin(pkin(5));
t137 = qJD(1) * t120;
t107 = (-t125 * t126 + t129 * t130) * t137;
t108 = (t125 * t130 + t126 * t129) * t137;
t124 = sin(qJ(4));
t128 = cos(qJ(4));
t100 = t107 * t124 + t108 * t128;
t138 = pkin(11) * t100;
t118 = cos(pkin(5)) * qJD(1);
t117 = t118 + qJD(2);
t113 = qJD(3) + t117;
t136 = pkin(1) * t118;
t116 = t130 * t136;
t135 = t126 * t137;
t104 = pkin(2) * t117 + t116 + (-pkin(8) - pkin(9)) * t135;
t134 = t130 * t137;
t110 = pkin(8) * t134 + t126 * t136;
t106 = pkin(9) * t134 + t110;
t97 = t129 * t104 - t106 * t125;
t93 = pkin(3) * t113 - pkin(10) * t108 + t97;
t98 = t125 * t104 + t129 * t106;
t95 = pkin(10) * t107 + t98;
t87 = t124 * t93 + t128 * t95;
t86 = -t124 * t95 + t128 * t93;
t119 = sin(pkin(6));
t121 = cos(pkin(6));
t112 = qJD(4) + t113;
t85 = pkin(4) * t112 - t121 * t138 + t86;
t111 = (-pkin(2) * t130 - pkin(1)) * t137;
t103 = -pkin(3) * t107 + t111;
t99 = t107 * t128 - t108 * t124;
t88 = -pkin(4) * t99 - t119 * t138 + t103;
t133 = t119 * t88 + t121 * t85;
t132 = t112 * t119 + t121 * t99;
t131 = qJD(1) ^ 2;
t127 = cos(qJ(5));
t123 = sin(qJ(5));
t109 = -pkin(8) * t135 + t116;
t96 = t112 * t121 - t119 * t99 + qJD(5);
t90 = t100 * t127 + t132 * t123;
t89 = -t100 * t123 + t132 * t127;
t84 = t132 * pkin(11) + t87;
t83 = -t119 * t85 + t121 * t88;
t82 = t133 * t123 + t127 * t84;
t81 = -t123 * t84 + t133 * t127;
t1 = t131 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t120 ^ 2 * t131 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + (-t103 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,2) * t99 / 0.2e1) * t99 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,3) * t96 / 0.2e1) * t96 + (t97 * mrSges(4,1) - t98 * mrSges(4,2) + Ifges(4,3) * t113 / 0.2e1) * t113 + (t83 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,5) * t96 + Ifges(6,1) * t90 / 0.2e1) * t90 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,6) * t99 + Ifges(5,3) * t112 / 0.2e1) * t112 + (t111 * mrSges(4,2) - t97 * mrSges(4,3) + Ifges(4,5) * t113 + Ifges(4,1) * t108 / 0.2e1) * t108 + (-t83 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,6) * t96 + Ifges(6,2) * t89 / 0.2e1) * t89 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t130 / 0.2e1) * t130 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t130 + Ifges(3,1) * t126 / 0.2e1) * t126) * t137 + (-t109 * t126 + t110 * t130) * mrSges(3,3)) * t137 + (t109 * mrSges(3,1) - t110 * mrSges(3,2) + Ifges(3,3) * t117 / 0.2e1 + (Ifges(3,5) * t126 + Ifges(3,6) * t130) * t137) * t117 + (-t111 * mrSges(4,1) + t98 * mrSges(4,3) + Ifges(4,4) * t108 + Ifges(4,6) * t113 + Ifges(4,2) * t107 / 0.2e1) * t107 + (t103 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,4) * t99 + Ifges(5,5) * t112 + Ifges(5,1) * t100 / 0.2e1) * t100;
T = t1;
