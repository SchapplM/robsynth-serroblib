% Calculate kinetic energy for
% S6RRPRPR5
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:42
% EndTime: 2019-03-09 10:28:43
% DurationCPUTime: 0.68s
% Computational Cost: add. (1164->115), mult. (3112->181), div. (0->0), fcn. (2514->12), ass. (0->52)
t133 = sin(pkin(11));
t136 = cos(pkin(11));
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t134 = sin(pkin(6));
t149 = qJD(1) * t134;
t122 = (t133 * t140 - t136 * t142) * t149;
t150 = cos(qJ(4));
t132 = sin(pkin(12));
t135 = cos(pkin(12));
t148 = cos(pkin(6)) * qJD(1);
t147 = pkin(1) * t148;
t130 = t142 * t147;
t131 = qJD(2) + t148;
t146 = t140 * t149;
t117 = pkin(2) * t131 + t130 + (-pkin(8) - qJ(3)) * t146;
t145 = t142 * t149;
t125 = pkin(8) * t145 + t140 * t147;
t120 = qJ(3) * t145 + t125;
t108 = t133 * t117 + t136 * t120;
t106 = pkin(9) * t131 + t108;
t123 = (t133 * t142 + t136 * t140) * t149;
t126 = qJD(3) + (-pkin(2) * t142 - pkin(1)) * t149;
t112 = pkin(3) * t122 - pkin(9) * t123 + t126;
t139 = sin(qJ(4));
t100 = t150 * t106 + t139 * t112;
t121 = qJD(4) + t122;
t95 = qJ(5) * t121 + t100;
t107 = t117 * t136 - t133 * t120;
t105 = -pkin(3) * t131 - t107;
t115 = t123 * t139 - t150 * t131;
t116 = t150 * t123 + t139 * t131;
t98 = pkin(4) * t115 - qJ(5) * t116 + t105;
t91 = t132 * t98 + t135 * t95;
t90 = -t132 * t95 + t135 * t98;
t99 = -t139 * t106 + t150 * t112;
t94 = -t121 * pkin(4) + qJD(5) - t99;
t143 = qJD(1) ^ 2;
t141 = cos(qJ(6));
t138 = sin(qJ(6));
t124 = -pkin(8) * t146 + t130;
t114 = qJD(6) + t115;
t110 = t116 * t135 + t121 * t132;
t109 = -t116 * t132 + t121 * t135;
t102 = t109 * t138 + t110 * t141;
t101 = t109 * t141 - t110 * t138;
t92 = -t109 * pkin(5) + t94;
t89 = pkin(10) * t109 + t91;
t88 = pkin(5) * t115 - pkin(10) * t110 + t90;
t87 = t138 * t88 + t141 * t89;
t86 = -t138 * t89 + t141 * t88;
t1 = m(3) * (pkin(1) ^ 2 * t134 ^ 2 * t143 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + t143 * Ifges(2,3) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t126 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t94 ^ 2) / 0.2e1 + (t126 * mrSges(4,2) - t107 * mrSges(4,3) + Ifges(4,1) * t123 / 0.2e1) * t123 + (t99 * mrSges(5,1) - t100 * mrSges(5,2) + Ifges(5,3) * t121 / 0.2e1) * t121 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,3) * t114 / 0.2e1) * t114 + (t94 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,1) * t110 / 0.2e1) * t110 - (-t126 * mrSges(4,1) + t108 * mrSges(4,3) + Ifges(4,4) * t123 - Ifges(4,2) * t122 / 0.2e1) * t122 + (t105 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,5) * t121 + Ifges(5,1) * t116 / 0.2e1) * t116 + (-t94 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t110 + Ifges(6,2) * t109 / 0.2e1) * t109 + (t92 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,5) * t114 + Ifges(7,1) * t102 / 0.2e1) * t102 + (-t92 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,4) * t102 + Ifges(7,6) * t114 + Ifges(7,2) * t101 / 0.2e1) * t101 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t142 / 0.2e1) * t142 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t142 + Ifges(3,1) * t140 / 0.2e1) * t140) * t149 + (-t124 * t140 + t125 * t142) * mrSges(3,3)) * t149 + (t124 * mrSges(3,1) + t107 * mrSges(4,1) - t125 * mrSges(3,2) - t108 * mrSges(4,2) + Ifges(4,5) * t123 - Ifges(4,6) * t122 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t131 + (Ifges(3,5) * t140 + Ifges(3,6) * t142) * t149) * t131 + (t105 * mrSges(5,1) + t90 * mrSges(6,1) - t91 * mrSges(6,2) - t100 * mrSges(5,3) - Ifges(5,4) * t116 + Ifges(6,5) * t110 - Ifges(5,6) * t121 + Ifges(6,6) * t109 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t115) * t115;
T  = t1;
