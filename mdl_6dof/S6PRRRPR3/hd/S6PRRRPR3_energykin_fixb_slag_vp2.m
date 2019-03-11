% Calculate kinetic energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:10:52
% EndTime: 2019-03-08 23:10:52
% DurationCPUTime: 0.39s
% Computational Cost: add. (382->96), mult. (817->143), div. (0->0), fcn. (553->10), ass. (0->42)
t129 = pkin(4) + pkin(10);
t128 = cos(qJ(4));
t115 = sin(qJ(4));
t117 = sin(qJ(2));
t112 = sin(pkin(6));
t127 = qJD(1) * t112;
t105 = qJD(2) * pkin(8) + t117 * t127;
t119 = cos(qJ(3));
t113 = cos(pkin(6));
t126 = qJD(1) * t113;
t108 = t119 * t126;
t116 = sin(qJ(3));
t94 = qJD(3) * pkin(3) + t108 + (-pkin(9) * qJD(2) - t105) * t116;
t125 = qJD(2) * t119;
t99 = t119 * t105 + t116 * t126;
t95 = pkin(9) * t125 + t99;
t89 = t115 * t94 + t128 * t95;
t120 = cos(qJ(2));
t124 = t120 * t127;
t111 = qJD(3) + qJD(4);
t86 = -qJ(5) * t111 - t89;
t88 = -t115 * t95 + t128 * t94;
t123 = qJD(5) - t88;
t103 = (t115 * t119 + t128 * t116) * qJD(2);
t100 = -t124 + (-pkin(3) * t119 - pkin(2)) * qJD(2);
t122 = -qJ(5) * t103 + t100;
t118 = cos(qJ(6));
t114 = sin(qJ(6));
t106 = -qJD(2) * pkin(2) - t124;
t102 = qJD(2) * t115 * t116 - t128 * t125;
t101 = qJD(6) + t103;
t98 = -t105 * t116 + t108;
t97 = t102 * t114 + t111 * t118;
t96 = t102 * t118 - t111 * t114;
t90 = pkin(4) * t102 + t122;
t87 = t129 * t102 + t122;
t85 = -t111 * pkin(4) + t123;
t84 = -pkin(5) * t102 - t86;
t83 = t103 * pkin(5) - t129 * t111 + t123;
t82 = t114 * t83 + t118 * t87;
t81 = -t114 * t87 + t118 * t83;
t1 = m(4) * (t106 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t90 ^ 2) / 0.2e1 + (t84 * mrSges(7,2) - t81 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (t98 * mrSges(4,1) - t99 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t84 * mrSges(7,1) + t82 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (m(3) * (t113 ^ 2 + (t117 ^ 2 + t120 ^ 2) * t112 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t81 * mrSges(7,1) - t82 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t101 / 0.2e1) * t101 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + t85 * mrSges(6,2) - t86 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t111) * t111 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t120 - mrSges(3,2) * t117) * t127 + (-t106 * mrSges(4,1) + t99 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t125 / 0.2e1) * t119 + (t106 * mrSges(4,2) - t98 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t119 + Ifges(4,1) * t116 / 0.2e1) * qJD(2)) * t116) * qJD(2) + (t85 * mrSges(6,1) + t100 * mrSges(5,2) - t88 * mrSges(5,3) - t90 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t103 + (-Ifges(6,4) + Ifges(5,5)) * t111) * t103 + (t100 * mrSges(5,1) + t86 * mrSges(6,1) - t90 * mrSges(6,2) - t89 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102 + (Ifges(6,5) - Ifges(5,6)) * t111 + (-Ifges(5,4) - Ifges(6,6)) * t103) * t102;
T  = t1;
