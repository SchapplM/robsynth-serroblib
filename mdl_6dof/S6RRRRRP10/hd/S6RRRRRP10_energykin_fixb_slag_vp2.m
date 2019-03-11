% Calculate kinetic energy for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:03
% EndTime: 2019-03-10 02:16:03
% DurationCPUTime: 0.56s
% Computational Cost: add. (1100->110), mult. (2404->166), div. (0->0), fcn. (1896->10), ass. (0->46)
t140 = cos(qJ(5));
t127 = sin(qJ(5));
t139 = cos(pkin(6)) * qJD(1);
t124 = qJD(2) + t139;
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t130 = sin(qJ(2));
t125 = sin(pkin(6));
t138 = t125 * qJD(1);
t136 = t130 * t138;
t116 = t124 * t129 + t132 * t136;
t133 = cos(qJ(2));
t135 = t133 * t138;
t120 = qJD(3) - t135;
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t106 = t116 * t131 + t120 * t128;
t115 = t124 * t132 - t129 * t136;
t114 = qJD(4) - t115;
t137 = pkin(1) * t139;
t118 = pkin(8) * t135 + t130 * t137;
t111 = pkin(9) * t124 + t118;
t113 = (-pkin(2) * t133 - pkin(9) * t130 - pkin(1)) * t138;
t103 = t132 * t111 + t129 * t113;
t101 = pkin(10) * t120 + t103;
t117 = -pkin(8) * t136 + t133 * t137;
t110 = -pkin(2) * t124 - t117;
t98 = -pkin(3) * t115 - pkin(10) * t116 + t110;
t91 = -t101 * t128 + t131 * t98;
t88 = pkin(4) * t114 - pkin(11) * t106 + t91;
t105 = -t116 * t128 + t120 * t131;
t92 = t131 * t101 + t128 * t98;
t90 = pkin(11) * t105 + t92;
t85 = t127 * t88 + t140 * t90;
t102 = -t129 * t111 + t113 * t132;
t100 = -pkin(3) * t120 - t102;
t84 = -t127 * t90 + t140 * t88;
t93 = -pkin(4) * t105 + t100;
t134 = qJD(1) ^ 2;
t112 = qJD(5) + t114;
t95 = t127 * t105 + t140 * t106;
t94 = -t140 * t105 + t106 * t127;
t86 = pkin(5) * t94 - qJ(6) * t95 + t93;
t83 = qJ(6) * t112 + t85;
t82 = -t112 * pkin(5) + qJD(6) - t84;
t1 = m(3) * (pkin(1) ^ 2 * t125 ^ 2 * t134 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + t134 * Ifges(2,3) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) / 0.2e1 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,3) * t120 / 0.2e1) * t120 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t114 / 0.2e1) * t114 + (t110 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,5) * t120 + Ifges(4,1) * t116 / 0.2e1) * t116 + (t100 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t114 + Ifges(5,1) * t106 / 0.2e1) * t106 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t133 / 0.2e1) * t133 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t133 + Ifges(3,1) * t130 / 0.2e1) * t130) * t138 + (-t117 * t130 + t118 * t133) * mrSges(3,3)) * t138 + (t117 * mrSges(3,1) - t118 * mrSges(3,2) + Ifges(3,3) * t124 / 0.2e1 + (Ifges(3,5) * t130 + Ifges(3,6) * t133) * t138) * t124 + (-t110 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t116 + Ifges(4,6) * t120 + Ifges(4,2) * t115 / 0.2e1) * t115 + (-t100 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,6) * t114 + Ifges(5,2) * t105 / 0.2e1) * t105 + (t93 * mrSges(6,2) + t82 * mrSges(7,2) - t84 * mrSges(6,3) - t86 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t95) * t95 + (t93 * mrSges(6,1) + t86 * mrSges(7,1) - t83 * mrSges(7,2) - t85 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t94 + (-Ifges(6,4) + Ifges(7,5)) * t95) * t94 + (t84 * mrSges(6,1) - t82 * mrSges(7,1) - t85 * mrSges(6,2) + t83 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t112 + (Ifges(7,4) + Ifges(6,5)) * t95 + (-Ifges(6,6) + Ifges(7,6)) * t94) * t112;
T  = t1;
