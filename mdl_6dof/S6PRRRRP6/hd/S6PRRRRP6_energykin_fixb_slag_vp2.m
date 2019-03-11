% Calculate kinetic energy for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:28
% EndTime: 2019-03-09 00:28:29
% DurationCPUTime: 0.48s
% Computational Cost: add. (620->97), mult. (1413->152), div. (0->0), fcn. (1105->12), ass. (0->44)
t138 = cos(qJ(2));
t129 = sin(pkin(6));
t145 = qJD(1) * t129;
t122 = qJD(2) * pkin(2) + t138 * t145;
t128 = sin(pkin(7));
t130 = cos(pkin(7));
t131 = cos(pkin(6));
t144 = qJD(1) * t131;
t148 = t122 * t130 + t128 * t144;
t135 = sin(qJ(2));
t143 = t128 * qJD(2);
t120 = pkin(9) * t143 + t135 * t145;
t134 = sin(qJ(3));
t137 = cos(qJ(3));
t107 = -t134 * t120 + t148 * t137;
t147 = cos(qJ(5));
t126 = qJD(2) * t130 + qJD(3);
t105 = -pkin(3) * t126 - t107;
t133 = sin(qJ(4));
t136 = cos(qJ(4));
t141 = t134 * t143;
t115 = t126 * t136 - t133 * t141;
t116 = t126 * t133 + t136 * t141;
t102 = -pkin(4) * t115 - pkin(11) * t116 + t105;
t132 = sin(qJ(5));
t108 = t137 * t120 + t148 * t134;
t106 = pkin(10) * t126 + t108;
t125 = t130 * t144;
t110 = t125 + (-t122 + (-pkin(3) * t137 - pkin(10) * t134) * qJD(2)) * t128;
t101 = t136 * t106 + t133 * t110;
t124 = -t137 * t143 + qJD(4);
t98 = pkin(11) * t124 + t101;
t94 = t132 * t102 + t147 * t98;
t100 = -t133 * t106 + t110 * t136;
t97 = -pkin(4) * t124 - t100;
t93 = t147 * t102 - t132 * t98;
t114 = qJD(5) - t115;
t113 = -t122 * t128 + t125;
t112 = t147 * t116 + t132 * t124;
t111 = t116 * t132 - t147 * t124;
t95 = pkin(5) * t111 - qJ(6) * t112 + t97;
t92 = qJ(6) * t114 + t94;
t91 = -t114 * pkin(5) + qJD(6) - t93;
t1 = m(4) * (t107 ^ 2 + t108 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t105 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t97 ^ 2) / 0.2e1 + (t107 * mrSges(4,1) - t108 * mrSges(4,2) + Ifges(4,3) * t126 / 0.2e1) * t126 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t124 / 0.2e1) * t124 + (t105 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,5) * t124 + Ifges(5,1) * t116 / 0.2e1) * t116 + (m(3) * (t131 ^ 2 + (t135 ^ 2 + t138 ^ 2) * t129 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t105 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t116 + Ifges(5,6) * t124 + Ifges(5,2) * t115 / 0.2e1) * t115 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t138 - mrSges(3,2) * t135) * t145 + (t113 * (-mrSges(4,1) * t137 + mrSges(4,2) * t134) + (Ifges(4,2) * t137 ^ 2 / 0.2e1 + (Ifges(4,4) * t137 + Ifges(4,1) * t134 / 0.2e1) * t134) * t143 + (-t107 * t134 + t108 * t137) * mrSges(4,3) + t126 * (Ifges(4,5) * t134 + Ifges(4,6) * t137)) * t128) * qJD(2) + (t93 * mrSges(6,1) - t91 * mrSges(7,1) - t94 * mrSges(6,2) + t92 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t114) * t114 + (t97 * mrSges(6,2) + t91 * mrSges(7,2) - t93 * mrSges(6,3) - t95 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t112 + (Ifges(7,4) + Ifges(6,5)) * t114) * t112 + (t97 * mrSges(6,1) + t95 * mrSges(7,1) - t92 * mrSges(7,2) - t94 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t111 + (-Ifges(6,6) + Ifges(7,6)) * t114 + (-Ifges(6,4) + Ifges(7,5)) * t112) * t111;
T  = t1;
