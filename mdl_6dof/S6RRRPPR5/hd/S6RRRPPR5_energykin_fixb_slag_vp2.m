% Calculate kinetic energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:11
% EndTime: 2019-03-09 15:38:11
% DurationCPUTime: 0.60s
% Computational Cost: add. (1374->114), mult. (3138->179), div. (0->0), fcn. (2538->12), ass. (0->51)
t143 = cos(qJ(2));
t140 = sin(qJ(2));
t135 = sin(pkin(6));
t149 = qJD(1) * t135;
t146 = t140 * t149;
t148 = cos(pkin(6)) * qJD(1);
t147 = pkin(1) * t148;
t126 = -pkin(8) * t146 + t143 * t147;
t132 = qJD(2) + t148;
t120 = -pkin(2) * t132 - t126;
t139 = sin(qJ(3));
t142 = cos(qJ(3));
t124 = t132 * t142 - t139 * t146;
t114 = -pkin(3) * t124 + qJD(4) + t120;
t125 = t132 * t139 + t142 * t146;
t134 = sin(pkin(11));
t150 = cos(pkin(11));
t115 = -t150 * t124 + t125 * t134;
t116 = t134 * t124 + t125 * t150;
t101 = pkin(4) * t115 - qJ(5) * t116 + t114;
t133 = sin(pkin(12));
t136 = cos(pkin(12));
t145 = t143 * t149;
t128 = qJD(3) - t145;
t127 = pkin(8) * t145 + t140 * t147;
t121 = pkin(9) * t132 + t127;
t123 = (-pkin(2) * t143 - pkin(9) * t140 - pkin(1)) * t149;
t111 = -t121 * t139 + t142 * t123;
t105 = pkin(3) * t128 - qJ(4) * t125 + t111;
t112 = t142 * t121 + t139 * t123;
t108 = qJ(4) * t124 + t112;
t98 = t134 * t105 + t150 * t108;
t96 = qJ(5) * t128 + t98;
t92 = t133 * t101 + t136 * t96;
t91 = t136 * t101 - t133 * t96;
t97 = t105 * t150 - t134 * t108;
t95 = -t128 * pkin(4) + qJD(5) - t97;
t144 = qJD(1) ^ 2;
t141 = cos(qJ(6));
t138 = sin(qJ(6));
t113 = qJD(6) + t115;
t110 = t116 * t136 + t128 * t133;
t109 = -t116 * t133 + t128 * t136;
t103 = t109 * t138 + t110 * t141;
t102 = t109 * t141 - t110 * t138;
t93 = -t109 * pkin(5) + t95;
t90 = pkin(10) * t109 + t92;
t89 = pkin(5) * t115 - pkin(10) * t110 + t91;
t88 = t138 * t89 + t141 * t90;
t87 = -t138 * t90 + t141 * t89;
t1 = t144 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t135 ^ 2 * t144 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + (t120 * mrSges(4,2) - t111 * mrSges(4,3) + Ifges(4,1) * t125 / 0.2e1) * t125 + (t114 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,1) * t116 / 0.2e1) * t116 + (t87 * mrSges(7,1) - t88 * mrSges(7,2) + Ifges(7,3) * t113 / 0.2e1) * t113 + (t95 * mrSges(6,2) - t91 * mrSges(6,3) + Ifges(6,1) * t110 / 0.2e1) * t110 + (-t120 * mrSges(4,1) + t112 * mrSges(4,3) + Ifges(4,4) * t125 + Ifges(4,2) * t124 / 0.2e1) * t124 + (-t95 * mrSges(6,1) + t92 * mrSges(6,3) + Ifges(6,4) * t110 + Ifges(6,2) * t109 / 0.2e1) * t109 + (t93 * mrSges(7,2) - t87 * mrSges(7,3) + Ifges(7,5) * t113 + Ifges(7,1) * t103 / 0.2e1) * t103 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t143 / 0.2e1) * t143 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t143 + Ifges(3,1) * t140 / 0.2e1) * t140) * t149 + (-t126 * t140 + t127 * t143) * mrSges(3,3)) * t149 + (t126 * mrSges(3,1) - t127 * mrSges(3,2) + Ifges(3,3) * t132 / 0.2e1 + (Ifges(3,5) * t140 + Ifges(3,6) * t143) * t149) * t132 + (-t93 * mrSges(7,1) + t88 * mrSges(7,3) + Ifges(7,4) * t103 + Ifges(7,6) * t113 + Ifges(7,2) * t102 / 0.2e1) * t102 + (t111 * mrSges(4,1) + t97 * mrSges(5,1) - t112 * mrSges(4,2) - t98 * mrSges(5,2) + Ifges(4,5) * t125 + Ifges(5,5) * t116 + Ifges(4,6) * t124 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t128) * t128 + (t114 * mrSges(5,1) + t91 * mrSges(6,1) - t92 * mrSges(6,2) - t98 * mrSges(5,3) - Ifges(5,4) * t116 + Ifges(6,5) * t110 - Ifges(5,6) * t128 + Ifges(6,6) * t109 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t115) * t115;
T  = t1;
