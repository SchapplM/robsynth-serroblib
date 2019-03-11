% Calculate kinetic energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:25
% EndTime: 2019-03-09 08:54:26
% DurationCPUTime: 0.58s
% Computational Cost: add. (1140->115), mult. (3100->180), div. (0->0), fcn. (2506->12), ass. (0->52)
t137 = sin(qJ(5));
t140 = cos(qJ(5));
t131 = sin(pkin(11));
t134 = cos(pkin(11));
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t132 = sin(pkin(6));
t147 = qJD(1) * t132;
t121 = (t131 * t141 + t134 * t138) * t147;
t146 = cos(pkin(6)) * qJD(1);
t129 = qJD(2) + t146;
t130 = sin(pkin(12));
t133 = cos(pkin(12));
t114 = t121 * t133 + t129 * t130;
t143 = t141 * t147;
t144 = t138 * t147;
t120 = t131 * t144 - t134 * t143;
t145 = pkin(1) * t146;
t128 = t141 * t145;
t115 = pkin(2) * t129 + t128 + (-pkin(8) - qJ(3)) * t144;
t123 = pkin(8) * t143 + t138 * t145;
t118 = qJ(3) * t143 + t123;
t108 = t131 * t115 + t134 * t118;
t106 = qJ(4) * t129 + t108;
t124 = qJD(3) + (-pkin(2) * t141 - pkin(1)) * t147;
t111 = pkin(3) * t120 - qJ(4) * t121 + t124;
t96 = -t106 * t130 + t133 * t111;
t93 = pkin(4) * t120 - pkin(9) * t114 + t96;
t113 = -t121 * t130 + t129 * t133;
t97 = t133 * t106 + t130 * t111;
t95 = pkin(9) * t113 + t97;
t90 = t137 * t93 + t140 * t95;
t107 = t115 * t134 - t131 * t118;
t89 = -t137 * t95 + t140 * t93;
t103 = t113 * t140 - t114 * t137;
t105 = -pkin(3) * t129 + qJD(4) - t107;
t98 = -pkin(4) * t113 + t105;
t142 = qJD(1) ^ 2;
t139 = cos(qJ(6));
t136 = sin(qJ(6));
t122 = -pkin(8) * t144 + t128;
t119 = qJD(5) + t120;
t104 = t113 * t137 + t114 * t140;
t101 = qJD(6) - t103;
t100 = t104 * t139 + t119 * t136;
t99 = -t104 * t136 + t119 * t139;
t91 = -pkin(5) * t103 - pkin(10) * t104 + t98;
t88 = pkin(10) * t119 + t90;
t87 = -pkin(5) * t119 - t89;
t86 = t136 * t91 + t139 * t88;
t85 = -t136 * t88 + t139 * t91;
t1 = m(5) * (t105 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t98 ^ 2) / 0.2e1 + t142 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t132 ^ 2 * t142 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + (-t87 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,2) * t99 / 0.2e1) * t99 + (t124 * mrSges(4,2) - t107 * mrSges(4,3) + Ifges(4,1) * t121 / 0.2e1) * t121 + (t89 * mrSges(6,1) - t90 * mrSges(6,2) + Ifges(6,3) * t119 / 0.2e1) * t119 + (t105 * mrSges(5,2) - t96 * mrSges(5,3) + Ifges(5,1) * t114 / 0.2e1) * t114 + (-t105 * mrSges(5,1) + t97 * mrSges(5,3) + Ifges(5,4) * t114 + Ifges(5,2) * t113 / 0.2e1) * t113 + (t98 * mrSges(6,2) - t89 * mrSges(6,3) + Ifges(6,5) * t119 + Ifges(6,1) * t104 / 0.2e1) * t104 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,6) * t99 + Ifges(7,3) * t101 / 0.2e1) * t101 + (-t98 * mrSges(6,1) + t90 * mrSges(6,3) + Ifges(6,4) * t104 + Ifges(6,6) * t119 + Ifges(6,2) * t103 / 0.2e1) * t103 + (t87 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,5) * t101 + Ifges(7,1) * t100 / 0.2e1) * t100 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t141 / 0.2e1) * t141 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t141 + Ifges(3,1) * t138 / 0.2e1) * t138) * t147 + (-t122 * t138 + t123 * t141) * mrSges(3,3)) * t147 + (t122 * mrSges(3,1) + t107 * mrSges(4,1) - t123 * mrSges(3,2) - t108 * mrSges(4,2) + Ifges(4,5) * t121 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t129 + (Ifges(3,5) * t138 + Ifges(3,6) * t141) * t147) * t129 + (t124 * mrSges(4,1) + t96 * mrSges(5,1) - t97 * mrSges(5,2) - t108 * mrSges(4,3) - Ifges(4,4) * t121 + Ifges(5,5) * t114 - Ifges(4,6) * t129 + Ifges(5,6) * t113 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t120) * t120;
T  = t1;
