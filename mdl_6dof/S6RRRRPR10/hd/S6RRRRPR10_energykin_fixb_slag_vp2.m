% Calculate kinetic energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:00:52
% EndTime: 2019-03-09 23:00:53
% DurationCPUTime: 0.57s
% Computational Cost: add. (962->111), mult. (2134->166), div. (0->0), fcn. (1666->10), ass. (0->49)
t140 = pkin(4) + pkin(11);
t139 = cos(qJ(4));
t125 = sin(qJ(4));
t137 = cos(pkin(6)) * qJD(1);
t121 = qJD(2) + t137;
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t127 = sin(qJ(2));
t122 = sin(pkin(6));
t138 = qJD(1) * t122;
t135 = t127 * t138;
t113 = t121 * t126 + t129 * t135;
t130 = cos(qJ(2));
t134 = t130 * t138;
t117 = qJD(3) - t134;
t136 = pkin(1) * t137;
t115 = pkin(8) * t134 + t127 * t136;
t109 = pkin(9) * t121 + t115;
t111 = (-pkin(2) * t130 - pkin(9) * t127 - pkin(1)) * t138;
t99 = -t109 * t126 + t111 * t129;
t93 = pkin(3) * t117 - pkin(10) * t113 + t99;
t100 = t109 * t129 + t111 * t126;
t112 = t121 * t129 - t126 * t135;
t96 = pkin(10) * t112 + t100;
t90 = t125 * t93 + t139 * t96;
t116 = qJD(4) + t117;
t88 = -qJ(5) * t116 - t90;
t89 = -t125 * t96 + t139 * t93;
t114 = -pkin(8) * t135 + t130 * t136;
t133 = qJD(5) - t89;
t103 = t112 * t125 + t113 * t139;
t108 = -pkin(2) * t121 - t114;
t104 = -pkin(3) * t112 + t108;
t132 = -qJ(5) * t103 + t104;
t131 = qJD(1) ^ 2;
t128 = cos(qJ(6));
t124 = sin(qJ(6));
t102 = -t112 * t139 + t113 * t125;
t101 = qJD(6) + t103;
t98 = t102 * t124 + t116 * t128;
t97 = t102 * t128 - t116 * t124;
t91 = pkin(4) * t102 + t132;
t87 = -pkin(4) * t116 + t133;
t86 = t102 * t140 + t132;
t85 = -pkin(5) * t102 - t88;
t84 = t103 * pkin(5) - t116 * t140 + t133;
t83 = t124 * t84 + t128 * t86;
t82 = -t124 * t86 + t128 * t84;
t1 = m(5) * (t104 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t87 ^ 2 + t88 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t85 ^ 2) / 0.2e1 + t131 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t122 ^ 2 * t131 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t108 ^ 2 + t99 ^ 2) / 0.2e1 + (t85 * mrSges(7,2) - t82 * mrSges(7,3) + Ifges(7,1) * t98 / 0.2e1) * t98 + (t99 * mrSges(4,1) - t100 * mrSges(4,2) + Ifges(4,3) * t117 / 0.2e1) * t117 + (-t85 * mrSges(7,1) + t83 * mrSges(7,3) + Ifges(7,4) * t98 + Ifges(7,2) * t97 / 0.2e1) * t97 + (t108 * mrSges(4,2) - t99 * mrSges(4,3) + Ifges(4,5) * t117 + Ifges(4,1) * t113 / 0.2e1) * t113 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t130 / 0.2e1) * t130 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t130 + Ifges(3,1) * t127 / 0.2e1) * t127) * t138 + (-t114 * t127 + t115 * t130) * mrSges(3,3)) * t138 + (t114 * mrSges(3,1) - t115 * mrSges(3,2) + Ifges(3,3) * t121 / 0.2e1 + (Ifges(3,5) * t127 + Ifges(3,6) * t130) * t138) * t121 + (-t108 * mrSges(4,1) + t100 * mrSges(4,3) + Ifges(4,4) * t113 + Ifges(4,6) * t117 + Ifges(4,2) * t112 / 0.2e1) * t112 + (t82 * mrSges(7,1) - t83 * mrSges(7,2) + Ifges(7,5) * t98 + Ifges(7,6) * t97 + Ifges(7,3) * t101 / 0.2e1) * t101 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + t87 * mrSges(6,2) - t88 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t116) * t116 + (t87 * mrSges(6,1) + t104 * mrSges(5,2) - t89 * mrSges(5,3) - t91 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t103 + (-Ifges(6,4) + Ifges(5,5)) * t116) * t103 + (t104 * mrSges(5,1) + t88 * mrSges(6,1) - t91 * mrSges(6,2) - t90 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102 + (Ifges(6,5) - Ifges(5,6)) * t116 + (-Ifges(5,4) - Ifges(6,6)) * t103) * t102;
T  = t1;
