% Calculate kinetic energy for
% S6RRRRRP9
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:12
% EndTime: 2019-03-10 02:00:12
% DurationCPUTime: 0.57s
% Computational Cost: add. (1110->110), mult. (2434->166), div. (0->0), fcn. (1926->10), ass. (0->46)
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t140 = cos(pkin(6)) * qJD(1);
t124 = qJD(2) + t140;
t129 = sin(qJ(3));
t133 = cos(qJ(3));
t130 = sin(qJ(2));
t125 = sin(pkin(6));
t139 = t125 * qJD(1);
t137 = t130 * t139;
t116 = t124 * t129 + t133 * t137;
t134 = cos(qJ(2));
t136 = t134 * t139;
t120 = qJD(3) - t136;
t128 = sin(qJ(4));
t132 = cos(qJ(4));
t106 = t116 * t132 + t120 * t128;
t115 = t124 * t133 - t129 * t137;
t114 = qJD(4) - t115;
t138 = pkin(1) * t140;
t118 = pkin(8) * t136 + t130 * t138;
t111 = pkin(9) * t124 + t118;
t113 = (-pkin(2) * t134 - pkin(9) * t130 - pkin(1)) * t139;
t104 = t133 * t111 + t129 * t113;
t102 = pkin(10) * t120 + t104;
t117 = -pkin(8) * t137 + t134 * t138;
t110 = -pkin(2) * t124 - t117;
t99 = -pkin(3) * t115 - pkin(10) * t116 + t110;
t92 = -t102 * t128 + t132 * t99;
t88 = pkin(4) * t114 - pkin(11) * t106 + t92;
t105 = -t116 * t128 + t120 * t132;
t93 = t132 * t102 + t128 * t99;
t91 = pkin(11) * t105 + t93;
t85 = t127 * t88 + t131 * t91;
t84 = -t127 * t91 + t131 * t88;
t103 = -t129 * t111 + t113 * t133;
t101 = -pkin(3) * t120 - t103;
t94 = -pkin(4) * t105 + t101;
t135 = qJD(1) ^ 2;
t112 = qJD(5) + t114;
t96 = t105 * t127 + t106 * t131;
t95 = t105 * t131 - t106 * t127;
t89 = -pkin(5) * t95 + qJD(6) + t94;
t83 = qJ(6) * t95 + t85;
t82 = pkin(5) * t112 - qJ(6) * t96 + t84;
t1 = m(6) * (t84 ^ 2 + t85 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t125 ^ 2 * t135 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + t135 * Ifges(2,3) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,3) * t120 / 0.2e1) * t120 + (t92 * mrSges(5,1) - t93 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (t110 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,5) * t120 + Ifges(4,1) * t116 / 0.2e1) * t116 + (t101 * mrSges(5,2) - t92 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t106 / 0.2e1) * t106 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t134 / 0.2e1) * t134 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t134 + Ifges(3,1) * t130 / 0.2e1) * t130) * t139 + (-t117 * t130 + t118 * t134) * mrSges(3,3)) * t139 + (t117 * mrSges(3,1) - t118 * mrSges(3,2) + Ifges(3,3) * t124 / 0.2e1 + (Ifges(3,5) * t130 + Ifges(3,6) * t134) * t139) * t124 + (-t110 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t116 + Ifges(4,6) * t120 + Ifges(4,2) * t115 / 0.2e1) * t115 + (-t101 * mrSges(5,1) + t93 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,6) * t114 + Ifges(5,2) * t105 / 0.2e1) * t105 + (t94 * mrSges(6,2) + t89 * mrSges(7,2) - t84 * mrSges(6,3) - t82 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t96) * t96 + (-t94 * mrSges(6,1) - t89 * mrSges(7,1) + t85 * mrSges(6,3) + t83 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t95 + (Ifges(6,4) + Ifges(7,4)) * t96) * t95 + (t84 * mrSges(6,1) + t82 * mrSges(7,1) - t85 * mrSges(6,2) - t83 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t112 + (Ifges(6,5) + Ifges(7,5)) * t96 + (Ifges(6,6) + Ifges(7,6)) * t95) * t112;
T  = t1;
