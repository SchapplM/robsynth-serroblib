% Calculate kinetic energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:32
% EndTime: 2019-03-09 16:11:33
% DurationCPUTime: 0.62s
% Computational Cost: add. (896->111), mult. (2012->164), div. (0->0), fcn. (1554->10), ass. (0->48)
t146 = -pkin(4) - pkin(5);
t145 = cos(qJ(3));
t132 = sin(qJ(2));
t134 = cos(qJ(2));
t128 = sin(pkin(6));
t142 = t128 * qJD(1);
t139 = t134 * t142;
t143 = cos(pkin(6)) * qJD(1);
t141 = pkin(1) * t143;
t119 = pkin(8) * t139 + t132 * t141;
t126 = qJD(2) + t143;
t112 = pkin(9) * t126 + t119;
t113 = (-pkin(2) * t134 - pkin(9) * t132 - pkin(1)) * t142;
t131 = sin(qJ(3));
t104 = t145 * t112 + t131 * t113;
t122 = qJD(3) - t139;
t102 = qJ(4) * t122 + t104;
t127 = sin(pkin(11));
t144 = cos(pkin(11));
t140 = t132 * t142;
t118 = -pkin(8) * t140 + t134 * t141;
t111 = -pkin(2) * t126 - t118;
t116 = -t145 * t126 + t131 * t140;
t117 = t131 * t126 + t145 * t140;
t98 = pkin(3) * t116 - qJ(4) * t117 + t111;
t94 = t144 * t102 + t127 * t98;
t103 = -t131 * t112 + t145 * t113;
t91 = t116 * qJ(5) + t94;
t93 = -t127 * t102 + t144 * t98;
t138 = pkin(3) * t122 - qJD(4) + t103;
t137 = qJD(5) - t93;
t106 = t144 * t117 + t127 * t122;
t136 = qJ(5) * t106 + t138;
t135 = qJD(1) ^ 2;
t133 = cos(qJ(6));
t130 = sin(qJ(6));
t115 = qJD(6) - t116;
t105 = t117 * t127 - t144 * t122;
t96 = t105 * t130 + t106 * t133;
t95 = t105 * t133 - t106 * t130;
t92 = pkin(4) * t105 - t136;
t90 = -t116 * pkin(4) + t137;
t89 = t146 * t105 + t136;
t88 = pkin(10) * t105 + t91;
t87 = -t106 * pkin(10) + t146 * t116 + t137;
t86 = t130 * t87 + t133 * t88;
t85 = -t130 * t88 + t133 * t87;
t1 = m(5) * (t138 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t111 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t128 ^ 2 * t135 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t135 * Ifges(2,3) / 0.2e1 + (t89 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,1) * t96 / 0.2e1) * t96 + (t103 * mrSges(4,1) - t104 * mrSges(4,2) + Ifges(4,3) * t122 / 0.2e1) * t122 + (-t89 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,4) * t96 + Ifges(7,2) * t95 / 0.2e1) * t95 + (t111 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,5) * t122 + Ifges(4,1) * t117 / 0.2e1) * t117 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t134 / 0.2e1) * t134 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t134 + Ifges(3,1) * t132 / 0.2e1) * t132) * t142 + (-t118 * t132 + t119 * t134) * mrSges(3,3)) * t142 + (t118 * mrSges(3,1) - t119 * mrSges(3,2) + Ifges(3,3) * t126 / 0.2e1 + (Ifges(3,5) * t132 + Ifges(3,6) * t134) * t142) * t126 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,5) * t96 + Ifges(7,6) * t95 + Ifges(7,3) * t115 / 0.2e1) * t115 + (-t138 * mrSges(5,2) + t90 * mrSges(6,2) - t93 * mrSges(5,3) - t92 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t106) * t106 + (-t138 * mrSges(5,1) + t92 * mrSges(6,1) - t91 * mrSges(6,2) - t94 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t105 + (-Ifges(5,4) + Ifges(6,5)) * t106) * t105 + (t111 * mrSges(4,1) + t93 * mrSges(5,1) - t90 * mrSges(6,1) - t94 * mrSges(5,2) - t104 * mrSges(4,3) + t91 * mrSges(6,3) - Ifges(4,4) * t117 - Ifges(4,6) * t122 + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t116 + (Ifges(6,4) + Ifges(5,5)) * t106 + (-Ifges(5,6) + Ifges(6,6)) * t105) * t116;
T  = t1;
