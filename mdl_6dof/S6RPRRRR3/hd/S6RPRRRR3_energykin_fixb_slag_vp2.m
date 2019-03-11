% Calculate kinetic energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:13
% EndTime: 2019-03-09 07:00:13
% DurationCPUTime: 0.58s
% Computational Cost: add. (618->97), mult. (1269->156), div. (0->0), fcn. (850->10), ass. (0->42)
t127 = m(3) / 0.2e1;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t118 = sin(qJ(3));
t126 = qJD(1) * t118;
t106 = qJD(3) * t117 + t121 * t126;
t122 = cos(qJ(3));
t111 = -qJD(1) * t122 + qJD(4);
t114 = cos(pkin(11));
t125 = -pkin(1) * t114 - pkin(2);
t100 = (-pkin(3) * t122 - pkin(8) * t118 + t125) * qJD(1);
t113 = sin(pkin(11));
t108 = (pkin(1) * t113 + pkin(7)) * qJD(1);
t102 = t118 * qJD(2) + t122 * t108;
t99 = qJD(3) * pkin(8) + t102;
t90 = t121 * t100 - t117 * t99;
t86 = pkin(4) * t111 - pkin(9) * t106 + t90;
t105 = qJD(3) * t121 - t117 * t126;
t91 = t117 * t100 + t121 * t99;
t89 = pkin(9) * t105 + t91;
t81 = t116 * t86 + t120 * t89;
t80 = -t116 * t89 + t120 * t86;
t101 = qJD(2) * t122 - t118 * t108;
t110 = qJD(5) + t111;
t98 = -qJD(3) * pkin(3) - t101;
t92 = -pkin(4) * t105 + t98;
t119 = cos(qJ(6));
t115 = sin(qJ(6));
t109 = t125 * qJD(1);
t107 = qJD(6) + t110;
t94 = t105 * t116 + t106 * t120;
t93 = t105 * t120 - t106 * t116;
t87 = -pkin(5) * t93 + t92;
t85 = t115 * t93 + t119 * t94;
t84 = -t115 * t94 + t119 * t93;
t79 = pkin(10) * t93 + t81;
t78 = pkin(5) * t110 - pkin(10) * t94 + t80;
t77 = t115 * t78 + t119 * t79;
t76 = -t115 * t79 + t119 * t78;
t1 = m(7) * (t76 ^ 2 + t77 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t109 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t127 + (t92 * mrSges(6,2) - t80 * mrSges(6,3) + Ifges(6,1) * t94 / 0.2e1) * t94 + (t87 * mrSges(7,2) - t76 * mrSges(7,3) + Ifges(7,1) * t85 / 0.2e1) * t85 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * t111 / 0.2e1) * t111 + (-t92 * mrSges(6,1) + t81 * mrSges(6,3) + Ifges(6,4) * t94 + Ifges(6,2) * t93 / 0.2e1) * t93 + (-t87 * mrSges(7,1) + t77 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t98 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * t111 + Ifges(5,1) * t106 / 0.2e1) * t106 + (t80 * mrSges(6,1) - t81 * mrSges(6,2) + Ifges(6,5) * t94 + Ifges(6,6) * t93 + Ifges(6,3) * t110 / 0.2e1) * t110 + (t76 * mrSges(7,1) - t77 * mrSges(7,2) + Ifges(7,5) * t85 + Ifges(7,6) * t84 + Ifges(7,3) * t107 / 0.2e1) * t107 + (-t98 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t106 + Ifges(5,6) * t111 + Ifges(5,2) * t105 / 0.2e1) * t105 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t109 * (-mrSges(4,1) * t122 + mrSges(4,2) * t118) + (-t101 * t118 + t102 * t122) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t118 + Ifges(4,6) * t122) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t114 * mrSges(3,1) - t113 * mrSges(3,2) + (t113 ^ 2 + t114 ^ 2) * t127 * pkin(1)) * pkin(1) + Ifges(4,2) * t122 ^ 2 / 0.2e1 + (Ifges(4,4) * t122 + Ifges(4,1) * t118 / 0.2e1) * t118) * qJD(1)) * qJD(1);
T  = t1;
