% Calculate kinetic energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:24
% EndTime: 2019-03-10 01:34:24
% DurationCPUTime: 0.57s
% Computational Cost: add. (1102->110), mult. (2428->166), div. (0->0), fcn. (1922->10), ass. (0->46)
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t132 = cos(qJ(2));
t123 = sin(pkin(6));
t137 = t123 * qJD(1);
t134 = t132 * t137;
t118 = qJD(3) - t134;
t117 = qJD(4) + t118;
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t128 = sin(qJ(2));
t138 = cos(pkin(6)) * qJD(1);
t136 = pkin(1) * t138;
t116 = pkin(8) * t134 + t128 * t136;
t122 = qJD(2) + t138;
t111 = pkin(9) * t122 + t116;
t112 = (-pkin(2) * t132 - pkin(9) * t128 - pkin(1)) * t137;
t127 = sin(qJ(3));
t131 = cos(qJ(3));
t101 = -t111 * t127 + t131 * t112;
t135 = t128 * t137;
t114 = t122 * t127 + t131 * t135;
t95 = pkin(3) * t118 - pkin(10) * t114 + t101;
t102 = t131 * t111 + t127 * t112;
t113 = t122 * t131 - t127 * t135;
t98 = pkin(10) * t113 + t102;
t90 = t126 * t95 + t130 * t98;
t88 = pkin(11) * t117 + t90;
t104 = t113 * t130 - t114 * t126;
t105 = t113 * t126 + t114 * t130;
t115 = -pkin(8) * t135 + t132 * t136;
t110 = -pkin(2) * t122 - t115;
t106 = -pkin(3) * t113 + t110;
t93 = -pkin(4) * t104 - pkin(11) * t105 + t106;
t84 = t125 * t93 + t129 * t88;
t83 = -t125 * t88 + t129 * t93;
t89 = -t126 * t98 + t130 * t95;
t87 = -pkin(4) * t117 - t89;
t133 = qJD(1) ^ 2;
t103 = qJD(5) - t104;
t100 = t105 * t129 + t117 * t125;
t99 = -t105 * t125 + t117 * t129;
t85 = -pkin(5) * t99 + qJD(6) + t87;
t82 = qJ(6) * t99 + t84;
t81 = pkin(5) * t103 - qJ(6) * t100 + t83;
t1 = m(7) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + t133 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t123 ^ 2 * t133 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * t118 / 0.2e1) * t118 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + Ifges(5,3) * t117 / 0.2e1) * t117 + (t110 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,5) * t118 + Ifges(4,1) * t114 / 0.2e1) * t114 + (t106 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,5) * t117 + Ifges(5,1) * t105 / 0.2e1) * t105 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t132 / 0.2e1) * t132 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t132 + Ifges(3,1) * t128 / 0.2e1) * t128) * t137 + (-t115 * t128 + t116 * t132) * mrSges(3,3)) * t137 + (t115 * mrSges(3,1) - t116 * mrSges(3,2) + Ifges(3,3) * t122 / 0.2e1 + (Ifges(3,5) * t128 + Ifges(3,6) * t132) * t137) * t122 + (-t110 * mrSges(4,1) + t102 * mrSges(4,3) + Ifges(4,4) * t114 + Ifges(4,6) * t118 + Ifges(4,2) * t113 / 0.2e1) * t113 + (-t106 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,6) * t117 + Ifges(5,2) * t104 / 0.2e1) * t104 + (-t87 * mrSges(6,1) - t85 * mrSges(7,1) + t84 * mrSges(6,3) + t82 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t99) * t99 + (t83 * mrSges(6,1) + t81 * mrSges(7,1) - t84 * mrSges(6,2) - t82 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t103 + (Ifges(6,6) + Ifges(7,6)) * t99) * t103 + (t87 * mrSges(6,2) + t85 * mrSges(7,2) - t83 * mrSges(6,3) - t81 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t100 + (Ifges(6,4) + Ifges(7,4)) * t99 + (Ifges(6,5) + Ifges(7,5)) * t103) * t100;
T  = t1;
