% Calculate kinetic energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:37:45
% EndTime: 2019-03-09 02:37:46
% DurationCPUTime: 0.50s
% Computational Cost: add. (526->98), mult. (1185->154), div. (0->0), fcn. (798->10), ass. (0->41)
t124 = m(3) / 0.2e1;
t110 = sin(pkin(11));
t113 = cos(pkin(11));
t111 = sin(pkin(10));
t123 = cos(pkin(10));
t112 = sin(pkin(9));
t105 = (pkin(1) * t112 + pkin(7)) * qJD(1);
t118 = cos(qJ(3));
t109 = t118 * qJD(2);
t116 = sin(qJ(3));
t94 = qJD(3) * pkin(3) + t109 + (-qJ(4) * qJD(1) - t105) * t116;
t122 = qJD(1) * t118;
t99 = t116 * qJD(2) + t118 * t105;
t97 = qJ(4) * t122 + t99;
t85 = t111 * t94 + t123 * t97;
t83 = qJD(3) * qJ(5) + t85;
t114 = cos(pkin(9));
t121 = -pkin(1) * t114 - pkin(2);
t101 = qJD(4) + (-pkin(3) * t118 + t121) * qJD(1);
t102 = qJD(1) * t111 * t116 - t123 * t122;
t103 = (t111 * t118 + t123 * t116) * qJD(1);
t90 = pkin(4) * t102 - qJ(5) * t103 + t101;
t79 = t110 * t90 + t113 * t83;
t78 = -t110 * t83 + t113 * t90;
t84 = -t111 * t97 + t123 * t94;
t82 = -qJD(3) * pkin(4) + qJD(5) - t84;
t117 = cos(qJ(6));
t115 = sin(qJ(6));
t106 = t121 * qJD(1);
t100 = qJD(6) + t102;
t98 = -t105 * t116 + t109;
t96 = qJD(3) * t110 + t103 * t113;
t95 = qJD(3) * t113 - t103 * t110;
t87 = t115 * t95 + t117 * t96;
t86 = -t115 * t96 + t117 * t95;
t80 = -t95 * pkin(5) + t82;
t77 = pkin(8) * t95 + t79;
t76 = pkin(5) * t102 - pkin(8) * t96 + t78;
t75 = t115 * t76 + t117 * t77;
t74 = -t115 * t77 + t117 * t76;
t1 = m(7) * (t74 ^ 2 + t75 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t82 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t124 + m(4) * (t106 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + (t82 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t96 / 0.2e1) * t96 + (t80 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,1) * t87 / 0.2e1) * t87 + (t101 * mrSges(5,2) - t84 * mrSges(5,3) + Ifges(5,1) * t103 / 0.2e1) * t103 + (-t82 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t96 + Ifges(6,2) * t95 / 0.2e1) * t95 + (-t80 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t87 + Ifges(7,2) * t86 / 0.2e1) * t86 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t101 * mrSges(5,1) + t78 * mrSges(6,1) - t79 * mrSges(6,2) - t85 * mrSges(5,3) - Ifges(5,4) * t103 + Ifges(6,5) * t96 + Ifges(6,6) * t95 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t102) * t102 + (t98 * mrSges(4,1) + t84 * mrSges(5,1) - t99 * mrSges(4,2) - t85 * mrSges(5,2) + Ifges(5,5) * t103 - Ifges(5,6) * t102 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3)) * qJD(3) + (t106 * (-mrSges(4,1) * t118 + mrSges(4,2) * t116) + (-t116 * t98 + t118 * t99) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t116 + Ifges(4,6) * t118) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (mrSges(3,1) * t114 - mrSges(3,2) * t112 + (t112 ^ 2 + t114 ^ 2) * t124 * pkin(1)) * pkin(1) + Ifges(4,2) * t118 ^ 2 / 0.2e1 + (Ifges(4,4) * t118 + Ifges(4,1) * t116 / 0.2e1) * t116) * qJD(1)) * qJD(1);
T  = t1;
