% Calculate kinetic energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:57
% EndTime: 2019-03-08 18:55:57
% DurationCPUTime: 0.27s
% Computational Cost: add. (320->78), mult. (739->119), div. (0->0), fcn. (566->12), ass. (0->38)
t117 = cos(pkin(6)) * qJD(1) + qJD(2);
t122 = sin(pkin(7));
t125 = cos(pkin(7));
t124 = cos(pkin(12));
t123 = sin(pkin(6));
t139 = qJD(1) * t123;
t135 = t124 * t139;
t142 = t117 * t122 + t125 * t135;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t121 = sin(pkin(12));
t136 = t121 * t139;
t107 = -t128 * t136 + t130 * t142;
t141 = cos(qJ(5));
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t102 = (-pkin(4) * t129 - pkin(10) * t127 - pkin(3)) * qJD(3) - t107;
t126 = sin(qJ(5));
t108 = t142 * t128 + t130 * t136;
t106 = qJD(3) * pkin(9) + t108;
t110 = t125 * t117 - t122 * t135;
t100 = t129 * t106 + t127 * t110;
t98 = qJD(4) * pkin(10) + t100;
t94 = t126 * t102 + t141 * t98;
t138 = qJD(3) * t127;
t137 = t129 * qJD(3);
t99 = -t127 * t106 + t129 * t110;
t97 = -qJD(4) * pkin(4) - t99;
t93 = t102 * t141 - t126 * t98;
t131 = qJD(1) ^ 2;
t118 = qJD(5) - t137;
t114 = t126 * qJD(4) + t138 * t141;
t113 = -qJD(4) * t141 + t126 * t138;
t105 = -qJD(3) * pkin(3) - t107;
t95 = t113 * pkin(5) - t114 * qJ(6) + t97;
t92 = t118 * qJ(6) + t94;
t91 = -t118 * pkin(5) + qJD(6) - t93;
t1 = m(3) * (t117 ^ 2 + (t121 ^ 2 + t124 ^ 2) * t131 * t123 ^ 2) / 0.2e1 + m(2) * t131 / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + (t99 * mrSges(5,1) - t100 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t93 * mrSges(6,1) - t91 * mrSges(7,1) - t94 * mrSges(6,2) + t92 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t118) * t118 + (t107 * mrSges(4,1) - t108 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t105 * mrSges(5,1) + t100 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t137 / 0.2e1) * t129 + (t105 * mrSges(5,2) - t99 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t129 + Ifges(5,1) * t127 / 0.2e1) * qJD(3)) * t127) * qJD(3) + (t97 * mrSges(6,2) + t91 * mrSges(7,2) - t93 * mrSges(6,3) - t95 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t114 + (Ifges(7,4) + Ifges(6,5)) * t118) * t114 + (t97 * mrSges(6,1) + t95 * mrSges(7,1) - t92 * mrSges(7,2) - t94 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t113 + (-Ifges(6,6) + Ifges(7,6)) * t118 + (-Ifges(6,4) + Ifges(7,5)) * t114) * t113;
T  = t1;
