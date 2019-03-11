% Calculate kinetic energy for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:17
% EndTime: 2019-03-09 20:12:18
% DurationCPUTime: 0.56s
% Computational Cost: add. (878->111), mult. (1908->166), div. (0->0), fcn. (1454->10), ass. (0->49)
t142 = pkin(3) + pkin(10);
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t140 = cos(pkin(6)) * qJD(1);
t123 = qJD(2) + t140;
t128 = sin(qJ(3));
t132 = cos(qJ(3));
t129 = sin(qJ(2));
t124 = sin(pkin(6));
t141 = qJD(1) * t124;
t138 = t129 * t141;
t114 = t123 * t128 + t132 * t138;
t133 = cos(qJ(2));
t137 = t133 * t141;
t119 = qJD(3) - t137;
t139 = pkin(1) * t140;
t116 = pkin(8) * t137 + t129 * t139;
t109 = pkin(9) * t123 + t116;
t111 = (-pkin(2) * t133 - pkin(9) * t129 - pkin(1)) * t141;
t101 = -t128 * t109 + t111 * t132;
t136 = qJD(4) - t101;
t92 = pkin(4) * t114 - t142 * t119 + t136;
t113 = -t132 * t123 + t128 * t138;
t115 = -pkin(8) * t138 + t133 * t139;
t108 = -pkin(2) * t123 - t115;
t135 = -qJ(4) * t114 + t108;
t94 = t142 * t113 + t135;
t88 = t127 * t92 + t131 * t94;
t102 = t132 * t109 + t128 * t111;
t100 = -t119 * qJ(4) - t102;
t87 = -t127 * t94 + t131 * t92;
t95 = -pkin(4) * t113 - t100;
t112 = qJD(5) + t114;
t134 = qJD(1) ^ 2;
t130 = cos(qJ(6));
t126 = sin(qJ(6));
t110 = qJD(6) + t112;
t104 = t113 * t127 + t119 * t131;
t103 = t113 * t131 - t119 * t127;
t99 = -pkin(3) * t119 + t136;
t98 = pkin(3) * t113 + t135;
t97 = t103 * t126 + t104 * t130;
t96 = t103 * t130 - t104 * t126;
t89 = -pkin(5) * t103 + t95;
t86 = pkin(11) * t103 + t88;
t85 = pkin(5) * t112 - pkin(11) * t104 + t87;
t84 = t126 * t85 + t130 * t86;
t83 = -t126 * t86 + t130 * t85;
t1 = m(3) * (pkin(1) ^ 2 * t124 ^ 2 * t134 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + t134 * Ifges(2,3) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t87 ^ 2 + t88 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t83 ^ 2 + t84 ^ 2 + t89 ^ 2) / 0.2e1 + (t89 * mrSges(7,2) - t83 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t87 * mrSges(6,1) - t88 * mrSges(6,2) + Ifges(6,3) * t112 / 0.2e1) * t112 + (-t89 * mrSges(7,1) + t84 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (t95 * mrSges(6,2) - t87 * mrSges(6,3) + Ifges(6,5) * t112 + Ifges(6,1) * t104 / 0.2e1) * t104 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t133 / 0.2e1) * t133 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t133 + Ifges(3,1) * t129 / 0.2e1) * t129) * t141 + (-t115 * t129 + t116 * t133) * mrSges(3,3)) * t141 + (t115 * mrSges(3,1) - t116 * mrSges(3,2) + Ifges(3,3) * t123 / 0.2e1 + (Ifges(3,5) * t129 + Ifges(3,6) * t133) * t141) * t123 + (t83 * mrSges(7,1) - t84 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t110 / 0.2e1) * t110 + (-t95 * mrSges(6,1) + t88 * mrSges(6,3) + Ifges(6,4) * t104 + Ifges(6,6) * t112 + Ifges(6,2) * t103 / 0.2e1) * t103 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + t99 * mrSges(5,2) - t100 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t119) * t119 + (t99 * mrSges(5,1) + t108 * mrSges(4,2) - t101 * mrSges(4,3) - t98 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t114 + (-Ifges(5,4) + Ifges(4,5)) * t119) * t114 + (t108 * mrSges(4,1) + t100 * mrSges(5,1) - t98 * mrSges(5,2) - t102 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t113 + (Ifges(5,5) - Ifges(4,6)) * t119 + (-Ifges(4,4) - Ifges(5,6)) * t114) * t113;
T  = t1;
