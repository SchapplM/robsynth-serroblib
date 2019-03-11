% Calculate kinetic energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:40
% EndTime: 2019-03-09 13:27:41
% DurationCPUTime: 0.63s
% Computational Cost: add. (1172->115), mult. (3100->182), div. (0->0), fcn. (2506->12), ass. (0->53)
t137 = sin(qJ(5));
t141 = cos(qJ(5));
t132 = sin(pkin(12));
t134 = cos(pkin(12));
t139 = sin(qJ(2));
t143 = cos(qJ(2));
t133 = sin(pkin(6));
t149 = qJD(1) * t133;
t123 = (t132 * t143 + t134 * t139) * t149;
t148 = cos(pkin(6)) * qJD(1);
t131 = qJD(2) + t148;
t138 = sin(qJ(4));
t142 = cos(qJ(4));
t115 = t123 * t142 + t131 * t138;
t145 = t143 * t149;
t146 = t139 * t149;
t122 = -t132 * t146 + t134 * t145;
t121 = qJD(4) - t122;
t147 = pkin(1) * t148;
t130 = t143 * t147;
t116 = pkin(2) * t131 + t130 + (-pkin(8) - qJ(3)) * t146;
t125 = pkin(8) * t145 + t139 * t147;
t119 = qJ(3) * t145 + t125;
t109 = t116 * t132 + t119 * t134;
t107 = pkin(9) * t131 + t109;
t126 = qJD(3) + (-pkin(2) * t143 - pkin(1)) * t149;
t112 = -pkin(3) * t122 - pkin(9) * t123 + t126;
t97 = -t107 * t138 + t112 * t142;
t94 = pkin(4) * t121 - pkin(10) * t115 + t97;
t114 = -t123 * t138 + t131 * t142;
t98 = t107 * t142 + t112 * t138;
t96 = pkin(10) * t114 + t98;
t91 = t137 * t94 + t141 * t96;
t108 = t116 * t134 - t119 * t132;
t90 = -t137 * t96 + t141 * t94;
t104 = t114 * t141 - t115 * t137;
t106 = -pkin(3) * t131 - t108;
t99 = -pkin(4) * t114 + t106;
t144 = qJD(1) ^ 2;
t140 = cos(qJ(6));
t136 = sin(qJ(6));
t124 = -pkin(8) * t146 + t130;
t120 = qJD(5) + t121;
t105 = t114 * t137 + t115 * t141;
t103 = qJD(6) - t104;
t101 = t105 * t140 + t120 * t136;
t100 = -t105 * t136 + t120 * t140;
t92 = -pkin(5) * t104 - pkin(11) * t105 + t99;
t89 = pkin(11) * t120 + t91;
t88 = -pkin(5) * t120 - t90;
t87 = t136 * t92 + t140 * t89;
t86 = -t136 * t89 + t140 * t92;
t1 = m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t126 ^ 2) / 0.2e1 + t144 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t133 ^ 2 * t144 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + (t126 * mrSges(4,2) - t108 * mrSges(4,3) + Ifges(4,1) * t123 / 0.2e1) * t123 + (t97 * mrSges(5,1) - t98 * mrSges(5,2) + Ifges(5,3) * t121 / 0.2e1) * t121 + (t90 * mrSges(6,1) - t91 * mrSges(6,2) + Ifges(6,3) * t120 / 0.2e1) * t120 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,3) * t103 / 0.2e1) * t103 + (-t126 * mrSges(4,1) + t109 * mrSges(4,3) + Ifges(4,4) * t123 + Ifges(4,2) * t122 / 0.2e1) * t122 + (t106 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,5) * t121 + Ifges(5,1) * t115 / 0.2e1) * t115 + (t99 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,5) * t120 + Ifges(6,1) * t105 / 0.2e1) * t105 + (t88 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,5) * t103 + Ifges(7,1) * t101 / 0.2e1) * t101 + (-t106 * mrSges(5,1) + t98 * mrSges(5,3) + Ifges(5,4) * t115 + Ifges(5,6) * t121 + Ifges(5,2) * t114 / 0.2e1) * t114 + (-t99 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,6) * t120 + Ifges(6,2) * t104 / 0.2e1) * t104 + (-t88 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,4) * t101 + Ifges(7,6) * t103 + Ifges(7,2) * t100 / 0.2e1) * t100 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t143 / 0.2e1) * t143 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t143 + Ifges(3,1) * t139 / 0.2e1) * t139) * t149 + (-t124 * t139 + t125 * t143) * mrSges(3,3)) * t149 + (t124 * mrSges(3,1) + t108 * mrSges(4,1) - t125 * mrSges(3,2) - t109 * mrSges(4,2) + Ifges(4,5) * t123 + Ifges(4,6) * t122 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t131 + (Ifges(3,5) * t139 + Ifges(3,6) * t143) * t149) * t131;
T  = t1;
