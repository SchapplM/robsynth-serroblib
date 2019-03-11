% Calculate kinetic energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:06
% EndTime: 2019-03-09 00:20:06
% DurationCPUTime: 0.49s
% Computational Cost: add. (622->97), mult. (1421->152), div. (0->0), fcn. (1113->12), ass. (0->44)
t139 = cos(qJ(2));
t129 = sin(pkin(6));
t146 = qJD(1) * t129;
t122 = qJD(2) * pkin(2) + t139 * t146;
t128 = sin(pkin(7));
t130 = cos(pkin(7));
t131 = cos(pkin(6));
t145 = qJD(1) * t131;
t148 = t122 * t130 + t128 * t145;
t135 = sin(qJ(2));
t144 = t128 * qJD(2);
t121 = pkin(9) * t144 + t135 * t146;
t134 = sin(qJ(3));
t138 = cos(qJ(3));
t108 = -t134 * t121 + t148 * t138;
t126 = qJD(2) * t130 + qJD(3);
t106 = -pkin(3) * t126 - t108;
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t142 = t134 * t144;
t116 = t126 * t137 - t133 * t142;
t117 = t126 * t133 + t137 * t142;
t103 = -pkin(4) * t116 - pkin(11) * t117 + t106;
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t109 = t138 * t121 + t148 * t134;
t107 = pkin(10) * t126 + t109;
t125 = t130 * t145;
t111 = t125 + (-t122 + (-pkin(3) * t138 - pkin(10) * t134) * qJD(2)) * t128;
t102 = t137 * t107 + t133 * t111;
t124 = -t138 * t144 + qJD(4);
t98 = pkin(11) * t124 + t102;
t94 = t132 * t103 + t136 * t98;
t93 = t136 * t103 - t132 * t98;
t101 = -t133 * t107 + t111 * t137;
t97 = -pkin(4) * t124 - t101;
t115 = qJD(5) - t116;
t114 = -t122 * t128 + t125;
t113 = t117 * t136 + t124 * t132;
t112 = -t117 * t132 + t124 * t136;
t95 = -pkin(5) * t112 + qJD(6) + t97;
t92 = qJ(6) * t112 + t94;
t91 = pkin(5) * t115 - qJ(6) * t113 + t93;
t1 = m(4) * (t108 ^ 2 + t109 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t102 ^ 2 + t106 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t97 ^ 2) / 0.2e1 + (t108 * mrSges(4,1) - t109 * mrSges(4,2) + Ifges(4,3) * t126 / 0.2e1) * t126 + (t101 * mrSges(5,1) - t102 * mrSges(5,2) + Ifges(5,3) * t124 / 0.2e1) * t124 + (t106 * mrSges(5,2) - t101 * mrSges(5,3) + Ifges(5,5) * t124 + Ifges(5,1) * t117 / 0.2e1) * t117 + (m(3) * (t131 ^ 2 + (t135 ^ 2 + t139 ^ 2) * t129 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t106 * mrSges(5,1) + t102 * mrSges(5,3) + Ifges(5,4) * t117 + Ifges(5,6) * t124 + Ifges(5,2) * t116 / 0.2e1) * t116 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t139 - mrSges(3,2) * t135) * t146 + (t114 * (-mrSges(4,1) * t138 + mrSges(4,2) * t134) + (Ifges(4,2) * t138 ^ 2 / 0.2e1 + (Ifges(4,4) * t138 + Ifges(4,1) * t134 / 0.2e1) * t134) * t144 + (-t108 * t134 + t109 * t138) * mrSges(4,3) + t126 * (Ifges(4,5) * t134 + Ifges(4,6) * t138)) * t128) * qJD(2) + (t93 * mrSges(6,1) + t91 * mrSges(7,1) - t94 * mrSges(6,2) - t92 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t115) * t115 + (t97 * mrSges(6,2) + t95 * mrSges(7,2) - t93 * mrSges(6,3) - t91 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t113 + (Ifges(6,5) + Ifges(7,5)) * t115) * t113 + (-t97 * mrSges(6,1) - t95 * mrSges(7,1) + t94 * mrSges(6,3) + t92 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t112 + (Ifges(6,6) + Ifges(7,6)) * t115 + (Ifges(6,4) + Ifges(7,4)) * t113) * t112;
T  = t1;
