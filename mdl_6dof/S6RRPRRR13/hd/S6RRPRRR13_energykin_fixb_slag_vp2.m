% Calculate kinetic energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:51
% EndTime: 2019-03-09 14:44:52
% DurationCPUTime: 0.52s
% Computational Cost: add. (814->113), mult. (1788->171), div. (0->0), fcn. (1320->10), ass. (0->50)
t142 = -pkin(2) - pkin(9);
t128 = sin(qJ(5));
t132 = cos(qJ(5));
t130 = sin(qJ(2));
t125 = sin(pkin(6));
t141 = qJD(1) * t125;
t137 = t130 * t141;
t118 = qJD(4) + t137;
t120 = pkin(8) * t137;
t126 = cos(pkin(6));
t140 = qJD(1) * t126;
t124 = qJD(2) + t140;
t134 = cos(qJ(2));
t100 = qJD(3) + t120 + t142 * t124 + (-pkin(1) * t126 * t134 + pkin(3) * t125 * t130) * qJD(1);
t136 = -qJ(3) * t130 - pkin(1);
t106 = (t142 * t134 + t136) * t141;
t129 = sin(qJ(4));
t133 = cos(qJ(4));
t93 = t129 * t100 + t133 * t106;
t91 = pkin(10) * t118 + t93;
t138 = t134 * t141;
t139 = pkin(1) * t140;
t115 = pkin(8) * t138 + t130 * t139;
t108 = -t124 * qJ(3) - t115;
t105 = pkin(3) * t138 - t108;
t112 = -t129 * t124 - t133 * t138;
t113 = t124 * t133 - t129 * t138;
t98 = -pkin(4) * t112 - pkin(10) * t113 + t105;
t87 = t128 * t98 + t132 * t91;
t111 = qJD(5) - t112;
t86 = -t128 * t91 + t132 * t98;
t92 = t100 * t133 - t129 * t106;
t114 = t134 * t139 - t120;
t90 = -pkin(4) * t118 - t92;
t135 = qJD(1) ^ 2;
t131 = cos(qJ(6));
t127 = sin(qJ(6));
t110 = (-pkin(2) * t134 + t136) * t141;
t109 = qJD(6) + t111;
t107 = -pkin(2) * t124 + qJD(3) - t114;
t102 = t113 * t132 + t118 * t128;
t101 = -t113 * t128 + t118 * t132;
t97 = t101 * t127 + t102 * t131;
t96 = t101 * t131 - t102 * t127;
t88 = -pkin(5) * t101 + t90;
t85 = pkin(11) * t101 + t87;
t84 = pkin(5) * t111 - pkin(11) * t102 + t86;
t83 = t127 * t84 + t131 * t85;
t82 = -t127 * t85 + t131 * t84;
t1 = m(7) * (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) / 0.2e1 + t135 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t125 ^ 2 * t135 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t110 ^ 2) / 0.2e1 + (t88 * mrSges(7,2) - t82 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t92 * mrSges(5,1) - t93 * mrSges(5,2) + Ifges(5,3) * t118 / 0.2e1) * t118 + (t86 * mrSges(6,1) - t87 * mrSges(6,2) + Ifges(6,3) * t111 / 0.2e1) * t111 + (-t88 * mrSges(7,1) + t83 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (t105 * mrSges(5,2) - t92 * mrSges(5,3) + Ifges(5,5) * t118 + Ifges(5,1) * t113 / 0.2e1) * t113 + (t90 * mrSges(6,2) - t86 * mrSges(6,3) + Ifges(6,5) * t111 + Ifges(6,1) * t102 / 0.2e1) * t102 + (-t105 * mrSges(5,1) + t93 * mrSges(5,3) + Ifges(5,4) * t113 + Ifges(5,6) * t118 + Ifges(5,2) * t112 / 0.2e1) * t112 + (t82 * mrSges(7,1) - t83 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t109 / 0.2e1) * t109 + (-t90 * mrSges(6,1) + t87 * mrSges(6,3) + Ifges(6,4) * t102 + Ifges(6,6) * t111 + Ifges(6,2) * t101 / 0.2e1) * t101 + ((-t108 * mrSges(4,1) + t110 * mrSges(4,2) + t115 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t134) * t141) * t134 + (t107 * mrSges(4,1) - t114 * mrSges(3,3) - t110 * mrSges(4,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t130 + (Ifges(3,4) + Ifges(4,6)) * t134) * t141) * t130) * t141 + (t114 * mrSges(3,1) - t115 * mrSges(3,2) + t107 * mrSges(4,2) - t108 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t124 + ((-Ifges(4,5) + Ifges(3,6)) * t134 + (-Ifges(4,4) + Ifges(3,5)) * t130) * t141) * t124;
T  = t1;
