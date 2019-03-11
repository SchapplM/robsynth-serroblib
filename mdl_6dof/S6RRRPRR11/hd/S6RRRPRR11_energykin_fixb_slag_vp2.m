% Calculate kinetic energy for
% S6RRRPRR11
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:44
% EndTime: 2019-03-09 19:24:45
% DurationCPUTime: 0.55s
% Computational Cost: add. (848->111), mult. (1860->166), div. (0->0), fcn. (1412->10), ass. (0->47)
t144 = cos(qJ(3));
t131 = sin(qJ(5));
t135 = cos(qJ(5));
t142 = cos(pkin(6)) * qJD(1);
t127 = qJD(2) + t142;
t132 = sin(qJ(3));
t133 = sin(qJ(2));
t128 = sin(pkin(6));
t143 = qJD(1) * t128;
t140 = t133 * t143;
t115 = t132 * t127 + t144 * t140;
t136 = cos(qJ(2));
t139 = t136 * t143;
t121 = qJD(3) - t139;
t141 = pkin(1) * t142;
t117 = pkin(8) * t139 + t133 * t141;
t110 = pkin(9) * t127 + t117;
t111 = (-pkin(2) * t136 - pkin(9) * t133 - pkin(1)) * t143;
t101 = -t132 * t110 + t144 * t111;
t138 = qJD(4) - t101;
t92 = -t115 * pkin(10) + (-pkin(3) - pkin(4)) * t121 + t138;
t114 = -t144 * t127 + t132 * t140;
t102 = t144 * t110 + t132 * t111;
t98 = t121 * qJ(4) + t102;
t95 = pkin(10) * t114 + t98;
t89 = t131 * t92 + t135 * t95;
t116 = -pkin(8) * t140 + t136 * t141;
t109 = -t127 * pkin(2) - t116;
t88 = -t131 * t95 + t135 * t92;
t96 = t114 * pkin(3) - t115 * qJ(4) + t109;
t104 = t114 * t135 - t115 * t131;
t93 = -pkin(4) * t114 - t96;
t137 = qJD(1) ^ 2;
t134 = cos(qJ(6));
t130 = sin(qJ(6));
t119 = qJD(5) - t121;
t105 = t114 * t131 + t115 * t135;
t103 = qJD(6) - t104;
t100 = t105 * t134 + t119 * t130;
t99 = -t105 * t130 + t119 * t134;
t97 = -t121 * pkin(3) + t138;
t90 = -pkin(5) * t104 - pkin(11) * t105 + t93;
t87 = pkin(11) * t119 + t89;
t86 = -pkin(5) * t119 - t88;
t85 = t130 * t90 + t134 * t87;
t84 = -t130 * t87 + t134 * t90;
t1 = m(4) * (t101 ^ 2 + t102 ^ 2 + t109 ^ 2) / 0.2e1 + t137 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t128 ^ 2 * t137 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + (-t86 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t88 * mrSges(6,1) - t89 * mrSges(6,2) + Ifges(6,3) * t119 / 0.2e1) * t119 + (t93 * mrSges(6,2) - t88 * mrSges(6,3) + Ifges(6,5) * t119 + Ifges(6,1) * t105 / 0.2e1) * t105 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t103 / 0.2e1) * t103 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t136 / 0.2e1) * t136 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t136 + Ifges(3,1) * t133 / 0.2e1) * t133) * t143 + (-t116 * t133 + t117 * t136) * mrSges(3,3)) * t143 + (t116 * mrSges(3,1) - t117 * mrSges(3,2) + Ifges(3,3) * t127 / 0.2e1 + (Ifges(3,5) * t133 + Ifges(3,6) * t136) * t143) * t127 + (-t93 * mrSges(6,1) + t89 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,6) * t119 + Ifges(6,2) * t104 / 0.2e1) * t104 + (t86 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t103 + Ifges(7,1) * t100 / 0.2e1) * t100 + (t101 * mrSges(4,1) - t97 * mrSges(5,1) - t102 * mrSges(4,2) + t98 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t121) * t121 + (t109 * mrSges(4,2) + t97 * mrSges(5,2) - t101 * mrSges(4,3) - t96 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t115 + (Ifges(5,4) + Ifges(4,5)) * t121) * t115 + (t109 * mrSges(4,1) + t96 * mrSges(5,1) - t98 * mrSges(5,2) - t102 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t114 + (-Ifges(4,6) + Ifges(5,6)) * t121 + (-Ifges(4,4) + Ifges(5,5)) * t115) * t114;
T  = t1;
