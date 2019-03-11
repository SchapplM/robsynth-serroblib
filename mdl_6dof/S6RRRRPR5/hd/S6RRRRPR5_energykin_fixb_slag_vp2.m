% Calculate kinetic energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:56
% EndTime: 2019-03-09 22:11:57
% DurationCPUTime: 0.65s
% Computational Cost: add. (682->107), mult. (1375->155), div. (0->0), fcn. (962->8), ass. (0->43)
t129 = -pkin(4) - pkin(5);
t128 = -pkin(8) - pkin(7);
t127 = pkin(7) * mrSges(3,3);
t126 = cos(qJ(4));
t113 = sin(qJ(4));
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t118 = cos(qJ(2));
t124 = t118 * qJD(1);
t115 = sin(qJ(2));
t125 = t115 * qJD(1);
t102 = -t114 * t125 + t117 * t124;
t103 = (t114 * t118 + t115 * t117) * qJD(1);
t108 = (-pkin(2) * t118 - pkin(1)) * qJD(1);
t88 = -pkin(3) * t102 - pkin(9) * t103 + t108;
t111 = qJD(2) + qJD(3);
t106 = qJD(2) * pkin(2) + t125 * t128;
t107 = t128 * t124;
t94 = t114 * t106 - t117 * t107;
t92 = pkin(9) * t111 + t94;
t84 = t113 * t88 + t126 * t92;
t93 = t117 * t106 + t114 * t107;
t101 = qJD(4) - t102;
t81 = t101 * qJ(5) + t84;
t123 = pkin(3) * t111 + t93;
t83 = -t113 * t92 + t126 * t88;
t122 = qJD(5) - t83;
t96 = t103 * t126 + t113 * t111;
t121 = qJ(5) * t96 + t123;
t116 = cos(qJ(6));
t112 = sin(qJ(6));
t98 = qJD(6) - t101;
t95 = t103 * t113 - t111 * t126;
t86 = t112 * t95 + t116 * t96;
t85 = -t112 * t96 + t116 * t95;
t82 = pkin(4) * t95 - t121;
t80 = -t101 * pkin(4) + t122;
t79 = t129 * t95 + t121;
t78 = pkin(10) * t95 + t81;
t77 = -t96 * pkin(10) + t101 * t129 + t122;
t76 = t112 * t77 + t116 * t78;
t75 = -t112 * t78 + t116 * t77;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t108 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t123 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,3) * t98 / 0.2e1) * t98 + (t93 * mrSges(4,1) - t94 * mrSges(4,2) + Ifges(4,3) * t111 / 0.2e1) * t111 + (t79 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,5) * t98 + Ifges(7,1) * t86 / 0.2e1) * t86 + (t108 * mrSges(4,2) - t93 * mrSges(4,3) + Ifges(4,5) * t111 + Ifges(4,1) * t103 / 0.2e1) * t103 + (-t79 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,6) * t98 + Ifges(7,2) * t85 / 0.2e1) * t85 + (-t108 * mrSges(4,1) + t94 * mrSges(4,3) + Ifges(4,4) * t103 + Ifges(4,6) * t111 + Ifges(4,2) * t102 / 0.2e1) * t102 + (-t123 * mrSges(5,2) + t80 * mrSges(6,2) - t83 * mrSges(5,3) - t82 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t96) * t96 + (-t123 * mrSges(5,1) + t82 * mrSges(6,1) - t81 * mrSges(6,2) - t84 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t95 + (-Ifges(5,4) + Ifges(6,5)) * t96) * t95 + (t83 * mrSges(5,1) - t80 * mrSges(6,1) - t84 * mrSges(5,2) + t81 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101 + (Ifges(6,4) + Ifges(5,5)) * t96 + (-Ifges(5,6) + Ifges(6,6)) * t95) * t101 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t115 ^ 2 + t118 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t127) * t118) * t118 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t118 + (Ifges(3,1) / 0.2e1 + t127) * t115) * t115) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t118 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t115) * qJD(2)) * qJD(1);
T  = t1;
