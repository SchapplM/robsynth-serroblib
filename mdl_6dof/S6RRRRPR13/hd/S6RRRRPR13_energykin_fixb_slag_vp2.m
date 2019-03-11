% Calculate kinetic energy for
% S6RRRRPR13
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:50:59
% EndTime: 2019-03-09 23:50:59
% DurationCPUTime: 0.68s
% Computational Cost: add. (924->111), mult. (2012->166), div. (0->0), fcn. (1554->10), ass. (0->49)
t148 = -pkin(4) - pkin(5);
t147 = cos(qJ(4));
t134 = sin(qJ(2));
t137 = cos(qJ(2));
t129 = sin(pkin(6));
t145 = t129 * qJD(1);
t142 = t137 * t145;
t146 = cos(pkin(6)) * qJD(1);
t144 = pkin(1) * t146;
t121 = pkin(8) * t142 + t134 * t144;
t128 = qJD(2) + t146;
t113 = pkin(9) * t128 + t121;
t116 = (-pkin(2) * t137 - pkin(9) * t134 - pkin(1)) * t145;
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t105 = t136 * t113 + t133 * t116;
t124 = qJD(3) - t142;
t103 = pkin(10) * t124 + t105;
t132 = sin(qJ(4));
t143 = t134 * t145;
t120 = -pkin(8) * t143 + t137 * t144;
t112 = -pkin(2) * t128 - t120;
t118 = t136 * t128 - t133 * t143;
t119 = t128 * t133 + t136 * t143;
t99 = -pkin(3) * t118 - pkin(10) * t119 + t112;
t95 = t147 * t103 + t132 * t99;
t104 = -t133 * t113 + t136 * t116;
t117 = qJD(4) - t118;
t92 = t117 * qJ(5) + t95;
t141 = pkin(3) * t124 + t104;
t94 = -t132 * t103 + t147 * t99;
t140 = qJD(5) - t94;
t107 = t119 * t147 + t132 * t124;
t139 = qJ(5) * t107 + t141;
t138 = qJD(1) ^ 2;
t135 = cos(qJ(6));
t131 = sin(qJ(6));
t115 = qJD(6) - t117;
t106 = t119 * t132 - t124 * t147;
t97 = t106 * t131 + t107 * t135;
t96 = t106 * t135 - t107 * t131;
t93 = pkin(4) * t106 - t139;
t91 = -t117 * pkin(4) + t140;
t90 = t106 * t148 + t139;
t89 = pkin(11) * t106 + t92;
t88 = -t107 * pkin(11) + t117 * t148 + t140;
t87 = t131 * t88 + t135 * t89;
t86 = -t131 * t89 + t135 * t88;
t1 = m(3) * (pkin(1) ^ 2 * t129 ^ 2 * t138 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + t138 * Ifges(2,3) / 0.2e1 + m(4) * (t104 ^ 2 + t105 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + (t90 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t104 * mrSges(4,1) - t105 * mrSges(4,2) + Ifges(4,3) * t124 / 0.2e1) * t124 + (-t90 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (t112 * mrSges(4,2) - t104 * mrSges(4,3) + Ifges(4,5) * t124 + Ifges(4,1) * t119 / 0.2e1) * t119 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t137 / 0.2e1) * t137 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t137 + Ifges(3,1) * t134 / 0.2e1) * t134) * t145 + (-t120 * t134 + t121 * t137) * mrSges(3,3)) * t145 + (t120 * mrSges(3,1) - t121 * mrSges(3,2) + Ifges(3,3) * t128 / 0.2e1 + (Ifges(3,5) * t134 + Ifges(3,6) * t137) * t145) * t128 + (-t112 * mrSges(4,1) + t105 * mrSges(4,3) + Ifges(4,4) * t119 + Ifges(4,6) * t124 + Ifges(4,2) * t118 / 0.2e1) * t118 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t115 / 0.2e1) * t115 + (t94 * mrSges(5,1) - t91 * mrSges(6,1) - t95 * mrSges(5,2) + t92 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t117) * t117 + (-t141 * mrSges(5,2) + t91 * mrSges(6,2) - t94 * mrSges(5,3) - t93 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t107 + (Ifges(6,4) + Ifges(5,5)) * t117) * t107 + (-t141 * mrSges(5,1) + t93 * mrSges(6,1) - t92 * mrSges(6,2) - t95 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t106 + (-Ifges(5,6) + Ifges(6,6)) * t117 + (-Ifges(5,4) + Ifges(6,5)) * t107) * t106;
T  = t1;
