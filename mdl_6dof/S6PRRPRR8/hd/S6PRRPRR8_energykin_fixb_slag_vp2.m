% Calculate kinetic energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:55
% EndTime: 2019-03-08 22:35:55
% DurationCPUTime: 0.49s
% Computational Cost: add. (482->103), mult. (1115->156), div. (0->0), fcn. (823->12), ass. (0->51)
t125 = sin(pkin(7));
t128 = cos(pkin(6));
t144 = qJD(1) * t128;
t136 = cos(qJ(2));
t126 = sin(pkin(6));
t145 = qJD(1) * t126;
t117 = qJD(2) * pkin(2) + t136 * t145;
t127 = cos(pkin(7));
t146 = t117 * t127;
t150 = t125 * t144 + t146;
t135 = cos(qJ(3));
t149 = t135 * t150;
t148 = -pkin(3) - pkin(10);
t122 = t127 * t144;
t131 = sin(qJ(3));
t147 = qJ(4) * t131;
t103 = t122 + (-t117 + (t135 * t148 - t147) * qJD(2)) * t125;
t130 = sin(qJ(5));
t134 = cos(qJ(5));
t123 = qJD(2) * t127 + qJD(3);
t132 = sin(qJ(2));
t143 = qJD(2) * t125;
t116 = pkin(9) * t143 + t132 * t145;
t114 = t131 * t116;
t141 = qJD(4) + t114;
t142 = qJD(2) * t131;
t98 = -t135 * t146 + (pkin(4) * t142 - t135 * t144) * t125 + t148 * t123 + t141;
t95 = t134 * t103 + t130 * t98;
t105 = t135 * t116 + t150 * t131;
t139 = t135 * t143;
t101 = -t123 * qJ(4) - t105;
t99 = pkin(4) * t139 - t101;
t94 = -t103 * t130 + t134 * t98;
t111 = -t123 * t130 - t134 * t139;
t133 = cos(qJ(6));
t129 = sin(qJ(6));
t119 = t125 * t142 + qJD(5);
t112 = t123 * t134 - t130 * t139;
t110 = qJD(6) - t111;
t109 = -t117 * t125 + t122;
t108 = t112 * t133 + t119 * t129;
t107 = -t112 * t129 + t119 * t133;
t106 = t122 + (-t117 + (-pkin(3) * t135 - t147) * qJD(2)) * t125;
t104 = -t114 + t149;
t100 = -pkin(3) * t123 + t141 - t149;
t96 = -pkin(5) * t111 - pkin(11) * t112 + t99;
t93 = pkin(11) * t119 + t95;
t92 = -pkin(5) * t119 - t94;
t91 = t129 * t96 + t133 * t93;
t90 = -t129 * t93 + t133 * t96;
t1 = m(5) * (t100 ^ 2 + t101 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t105 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t99 ^ 2) / 0.2e1 + (t94 * mrSges(6,1) - t95 * mrSges(6,2) + Ifges(6,3) * t119 / 0.2e1) * t119 + (t90 * mrSges(7,1) - t91 * mrSges(7,2) + Ifges(7,3) * t110 / 0.2e1) * t110 + (t99 * mrSges(6,2) - t94 * mrSges(6,3) + Ifges(6,5) * t119 + Ifges(6,1) * t112 / 0.2e1) * t112 + (t92 * mrSges(7,2) - t90 * mrSges(7,3) + Ifges(7,5) * t110 + Ifges(7,1) * t108 / 0.2e1) * t108 + (m(2) / 0.2e1 + m(3) * (t128 ^ 2 + (t132 ^ 2 + t136 ^ 2) * t126 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t99 * mrSges(6,1) + t95 * mrSges(6,3) + Ifges(6,4) * t112 + Ifges(6,6) * t119 + Ifges(6,2) * t111 / 0.2e1) * t111 + (-t92 * mrSges(7,1) + t91 * mrSges(7,3) + Ifges(7,4) * t108 + Ifges(7,6) * t110 + Ifges(7,2) * t107 / 0.2e1) * t107 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t136 - mrSges(3,2) * t132) * t145 + ((-t109 * mrSges(4,1) - t101 * mrSges(5,1) + t106 * mrSges(5,2) + t105 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t139) * t135 + (t100 * mrSges(5,1) + t109 * mrSges(4,2) - t104 * mrSges(4,3) - t106 * mrSges(5,3) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t131 + (Ifges(4,4) + Ifges(5,6)) * t135) * t143) * t131) * t125) * qJD(2) + (t104 * mrSges(4,1) - t105 * mrSges(4,2) + t100 * mrSges(5,2) - t101 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t123 + ((-Ifges(5,5) + Ifges(4,6)) * t135 + (-Ifges(5,4) + Ifges(4,5)) * t131) * t143) * t123;
T  = t1;
