% Calculate kinetic energy for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:52
% EndTime: 2019-03-09 16:53:53
% DurationCPUTime: 0.51s
% Computational Cost: add. (1064->110), mult. (2428->164), div. (0->0), fcn. (1922->10), ass. (0->45)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t130 = cos(qJ(2));
t122 = sin(pkin(6));
t136 = qJD(1) * t122;
t132 = t130 * t136;
t116 = qJD(3) - t132;
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t127 = sin(qJ(2));
t135 = cos(pkin(6)) * qJD(1);
t134 = pkin(1) * t135;
t115 = pkin(8) * t132 + t127 * t134;
t120 = qJD(2) + t135;
t110 = pkin(9) * t120 + t115;
t111 = (-pkin(2) * t130 - pkin(9) * t127 - pkin(1)) * t136;
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t100 = -t110 * t126 + t129 * t111;
t133 = t127 * t136;
t113 = t120 * t126 + t129 * t133;
t94 = pkin(3) * t116 - qJ(4) * t113 + t100;
t101 = t129 * t110 + t126 * t111;
t112 = t120 * t129 - t126 * t133;
t97 = qJ(4) * t112 + t101;
t89 = t121 * t94 + t123 * t97;
t87 = pkin(10) * t116 + t89;
t114 = -pkin(8) * t133 + t130 * t134;
t109 = -pkin(2) * t120 - t114;
t103 = -pkin(3) * t112 + qJD(4) + t109;
t104 = t112 * t123 - t113 * t121;
t105 = t112 * t121 + t113 * t123;
t92 = -pkin(4) * t104 - pkin(10) * t105 + t103;
t83 = t125 * t92 + t128 * t87;
t88 = -t121 * t97 + t123 * t94;
t82 = -t125 * t87 + t128 * t92;
t86 = -pkin(4) * t116 - t88;
t131 = qJD(1) ^ 2;
t102 = qJD(5) - t104;
t99 = t105 * t128 + t116 * t125;
t98 = -t105 * t125 + t116 * t128;
t84 = -pkin(5) * t98 + qJD(6) + t86;
t81 = qJ(6) * t98 + t83;
t80 = pkin(5) * t102 - qJ(6) * t99 + t82;
t1 = m(3) * (pkin(1) ^ 2 * t122 ^ 2 * t131 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t131 * Ifges(2,3) / 0.2e1 + m(5) * (t103 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) / 0.2e1 + (t109 * mrSges(4,2) - t100 * mrSges(4,3) + Ifges(4,1) * t113 / 0.2e1) * t113 + (t103 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (-t109 * mrSges(4,1) + t101 * mrSges(4,3) + Ifges(4,4) * t113 + Ifges(4,2) * t112 / 0.2e1) * t112 + (-t103 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t130 / 0.2e1) * t130 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t130 + Ifges(3,1) * t127 / 0.2e1) * t127) * t136 + (-t114 * t127 + t115 * t130) * mrSges(3,3)) * t136 + (t114 * mrSges(3,1) - t115 * mrSges(3,2) + Ifges(3,3) * t120 / 0.2e1 + (Ifges(3,5) * t127 + Ifges(3,6) * t130) * t136) * t120 + (t86 * mrSges(6,2) + t84 * mrSges(7,2) - t82 * mrSges(6,3) - t80 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t99) * t99 + (-t86 * mrSges(6,1) - t84 * mrSges(7,1) + t83 * mrSges(6,3) + t81 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t98 + (Ifges(6,4) + Ifges(7,4)) * t99) * t98 + (t100 * mrSges(4,1) + t88 * mrSges(5,1) - t101 * mrSges(4,2) - t89 * mrSges(5,2) + Ifges(4,5) * t113 + Ifges(5,5) * t105 + Ifges(4,6) * t112 + Ifges(5,6) * t104 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t116) * t116 + (t82 * mrSges(6,1) + t80 * mrSges(7,1) - t83 * mrSges(6,2) - t81 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t102 + (Ifges(6,5) + Ifges(7,5)) * t99 + (Ifges(6,6) + Ifges(7,6)) * t98) * t102;
T  = t1;
