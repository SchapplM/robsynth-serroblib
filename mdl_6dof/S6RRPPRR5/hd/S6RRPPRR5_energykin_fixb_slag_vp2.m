% Calculate kinetic energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:18
% EndTime: 2019-03-09 09:07:18
% DurationCPUTime: 0.42s
% Computational Cost: add. (480->106), mult. (1080->150), div. (0->0), fcn. (716->8), ass. (0->39)
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t122 = sin(qJ(2));
t125 = cos(qJ(2));
t118 = sin(pkin(6));
t131 = qJD(1) * t118;
t127 = t125 * t131;
t128 = t122 * t131;
t98 = -pkin(1) * t131 - pkin(2) * t127 - qJ(3) * t128;
t95 = pkin(3) * t127 + qJD(4) - t98;
t88 = (pkin(4) * t125 - pkin(9) * t122) * t131 + t95;
t130 = cos(pkin(6)) * qJD(1);
t116 = qJD(2) + t130;
t129 = pkin(1) * t130;
t103 = pkin(8) * t127 + t122 * t129;
t97 = t116 * qJ(3) + t103;
t94 = -qJ(4) * t127 + t97;
t91 = -pkin(9) * t116 + t94;
t84 = t121 * t88 + t124 * t91;
t102 = -pkin(8) * t128 + t125 * t129;
t96 = -t116 * pkin(2) + qJD(3) - t102;
t83 = -t121 * t91 + t124 * t88;
t90 = -t116 * pkin(3) - qJ(4) * t128 + t96;
t100 = -t116 * t124 - t121 * t128;
t87 = pkin(4) * t116 - t90;
t126 = qJD(1) ^ 2;
t123 = cos(qJ(6));
t120 = sin(qJ(6));
t104 = qJD(5) + t127;
t101 = -t116 * t121 + t124 * t128;
t99 = qJD(6) - t100;
t93 = t101 * t123 + t104 * t120;
t92 = -t101 * t120 + t104 * t123;
t85 = -pkin(5) * t100 - pkin(10) * t101 + t87;
t82 = pkin(10) * t104 + t84;
t81 = -pkin(5) * t104 - t83;
t80 = t120 * t85 + t123 * t82;
t79 = -t120 * t82 + t123 * t85;
t1 = m(3) * (pkin(1) ^ 2 * t118 ^ 2 * t126 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + Ifges(2,3) * t126 / 0.2e1 + m(5) * (t90 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + (t79 * mrSges(7,1) - t80 * mrSges(7,2) + Ifges(7,3) * t99 / 0.2e1) * t99 + (t83 * mrSges(6,1) - t84 * mrSges(6,2) + Ifges(6,3) * t104 / 0.2e1) * t104 + (t81 * mrSges(7,2) - t79 * mrSges(7,3) + Ifges(7,5) * t99 + Ifges(7,1) * t93 / 0.2e1) * t93 + (t87 * mrSges(6,2) - t83 * mrSges(6,3) + Ifges(6,5) * t104 + Ifges(6,1) * t101 / 0.2e1) * t101 + (-t81 * mrSges(7,1) + t80 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,6) * t99 + Ifges(7,2) * t92 / 0.2e1) * t92 + (-t87 * mrSges(6,1) + t84 * mrSges(6,3) + Ifges(6,4) * t101 + Ifges(6,6) * t104 + Ifges(6,2) * t100 / 0.2e1) * t100 + ((-t98 * mrSges(4,1) + t95 * mrSges(5,1) + t97 * mrSges(4,2) + t103 * mrSges(3,3) - t94 * mrSges(5,3) + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t125) * t131) * t125 + (t96 * mrSges(4,2) + t95 * mrSges(5,2) - t102 * mrSges(3,3) - t98 * mrSges(4,3) - t90 * mrSges(5,3) + (-pkin(1) * mrSges(3,2) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t122 + (Ifges(3,4) - Ifges(5,4) - Ifges(4,5)) * t125) * t131) * t122) * t131 + (t102 * mrSges(3,1) - t96 * mrSges(4,1) - t90 * mrSges(5,1) - t103 * mrSges(3,2) + t94 * mrSges(5,2) + t97 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t116 + ((Ifges(3,6) - Ifges(4,6) + Ifges(5,6)) * t125 + (Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * t122) * t131) * t116;
T  = t1;
