% Calculate kinetic energy for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:10
% EndTime: 2019-03-09 10:21:10
% DurationCPUTime: 0.69s
% Computational Cost: add. (1154->115), mult. (3100->181), div. (0->0), fcn. (2506->12), ass. (0->52)
t129 = sin(pkin(11));
t132 = cos(pkin(11));
t136 = sin(qJ(2));
t139 = cos(qJ(2));
t130 = sin(pkin(6));
t146 = qJD(1) * t130;
t119 = (t129 * t136 - t132 * t139) * t146;
t128 = sin(pkin(12));
t131 = cos(pkin(12));
t120 = (t129 * t139 + t132 * t136) * t146;
t145 = cos(pkin(6)) * qJD(1);
t127 = qJD(2) + t145;
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t113 = t120 * t138 + t127 * t135;
t118 = qJD(4) + t119;
t144 = pkin(1) * t145;
t126 = t139 * t144;
t143 = t136 * t146;
t114 = pkin(2) * t127 + t126 + (-pkin(8) - qJ(3)) * t143;
t142 = t139 * t146;
t122 = pkin(8) * t142 + t136 * t144;
t117 = qJ(3) * t142 + t122;
t107 = t129 * t114 + t132 * t117;
t105 = pkin(9) * t127 + t107;
t123 = qJD(3) + (-pkin(2) * t139 - pkin(1)) * t146;
t110 = pkin(3) * t119 - pkin(9) * t120 + t123;
t95 = -t105 * t135 + t138 * t110;
t92 = pkin(4) * t118 - qJ(5) * t113 + t95;
t112 = -t120 * t135 + t127 * t138;
t96 = t138 * t105 + t135 * t110;
t94 = qJ(5) * t112 + t96;
t89 = t128 * t92 + t131 * t94;
t106 = t114 * t132 - t129 * t117;
t88 = -t128 * t94 + t131 * t92;
t102 = t112 * t131 - t113 * t128;
t104 = -pkin(3) * t127 - t106;
t97 = -pkin(4) * t112 + qJD(5) + t104;
t140 = qJD(1) ^ 2;
t137 = cos(qJ(6));
t134 = sin(qJ(6));
t121 = -pkin(8) * t143 + t126;
t103 = t112 * t128 + t113 * t131;
t100 = qJD(6) - t102;
t99 = t103 * t137 + t118 * t134;
t98 = -t103 * t134 + t118 * t137;
t90 = -pkin(5) * t102 - pkin(10) * t103 + t97;
t87 = pkin(10) * t118 + t89;
t86 = -pkin(5) * t118 - t88;
t85 = t134 * t90 + t137 * t87;
t84 = -t134 * t87 + t137 * t90;
t1 = m(5) * (t104 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t130 ^ 2 * t140 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + t140 * Ifges(2,3) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t123 ^ 2) / 0.2e1 + (t86 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,1) * t99 / 0.2e1) * t99 + (t123 * mrSges(4,2) - t106 * mrSges(4,3) + Ifges(4,1) * t120 / 0.2e1) * t120 + (t104 * mrSges(5,2) - t95 * mrSges(5,3) + Ifges(5,1) * t113 / 0.2e1) * t113 + (t97 * mrSges(6,2) - t88 * mrSges(6,3) + Ifges(6,1) * t103 / 0.2e1) * t103 + (-t86 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,2) * t98 / 0.2e1) * t98 - (-t123 * mrSges(4,1) + t107 * mrSges(4,3) + Ifges(4,4) * t120 - Ifges(4,2) * t119 / 0.2e1) * t119 + (-t104 * mrSges(5,1) + t96 * mrSges(5,3) + Ifges(5,4) * t113 + Ifges(5,2) * t112 / 0.2e1) * t112 + (-t97 * mrSges(6,1) + t89 * mrSges(6,3) + Ifges(6,4) * t103 + Ifges(6,2) * t102 / 0.2e1) * t102 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,5) * t99 + Ifges(7,6) * t98 + Ifges(7,3) * t100 / 0.2e1) * t100 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t139 / 0.2e1) * t139 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t139 + Ifges(3,1) * t136 / 0.2e1) * t136) * t146 + (-t121 * t136 + t122 * t139) * mrSges(3,3)) * t146 + (t121 * mrSges(3,1) + t106 * mrSges(4,1) - t122 * mrSges(3,2) - t107 * mrSges(4,2) + Ifges(4,5) * t120 - Ifges(4,6) * t119 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t127 + (Ifges(3,5) * t136 + Ifges(3,6) * t139) * t146) * t127 + (t95 * mrSges(5,1) + t88 * mrSges(6,1) - t96 * mrSges(5,2) - t89 * mrSges(6,2) + Ifges(5,5) * t113 + Ifges(6,5) * t103 + Ifges(5,6) * t112 + Ifges(6,6) * t102 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t118) * t118;
T  = t1;
