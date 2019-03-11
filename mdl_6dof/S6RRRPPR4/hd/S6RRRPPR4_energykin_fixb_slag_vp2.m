% Calculate kinetic energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:50
% EndTime: 2019-03-09 15:31:51
% DurationCPUTime: 0.59s
% Computational Cost: add. (684->108), mult. (1445->154), div. (0->0), fcn. (1002->8), ass. (0->41)
t125 = -pkin(4) - pkin(5);
t124 = pkin(7) * mrSges(3,3);
t109 = sin(pkin(10));
t123 = cos(pkin(10));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t112 = sin(qJ(2));
t122 = t112 * qJD(1);
t101 = qJD(2) * t111 + t114 * t122;
t115 = cos(qJ(2));
t121 = t115 * qJD(1);
t107 = qJD(3) - t121;
t104 = pkin(7) * t121 + qJD(2) * pkin(8);
t99 = (-pkin(2) * t115 - pkin(8) * t112 - pkin(1)) * qJD(1);
t93 = -t104 * t111 + t114 * t99;
t87 = pkin(3) * t107 - qJ(4) * t101 + t93;
t100 = qJD(2) * t114 - t111 * t122;
t94 = t114 * t104 + t111 * t99;
t90 = qJ(4) * t100 + t94;
t82 = t109 * t87 + t123 * t90;
t80 = t107 * qJ(5) + t82;
t103 = -qJD(2) * pkin(2) + pkin(7) * t122;
t81 = -t109 * t90 + t123 * t87;
t120 = qJD(5) - t81;
t119 = pkin(3) * t100 - qJD(4) - t103;
t92 = t109 * t100 + t101 * t123;
t118 = qJ(5) * t92 + t119;
t113 = cos(qJ(6));
t110 = sin(qJ(6));
t106 = qJD(6) - t107;
t91 = -t100 * t123 + t101 * t109;
t85 = t110 * t91 + t113 * t92;
t84 = -t110 * t92 + t113 * t91;
t83 = pkin(4) * t91 - t118;
t79 = -t107 * pkin(4) + t120;
t78 = t125 * t91 + t118;
t77 = pkin(9) * t91 + t80;
t76 = -t92 * pkin(9) + t107 * t125 + t120;
t75 = t110 * t76 + t113 * t77;
t74 = -t110 * t77 + t113 * t76;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t103 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) / 0.2e1 + (t78 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,1) * t85 / 0.2e1) * t85 + (t103 * mrSges(4,2) - t93 * mrSges(4,3) + Ifges(4,1) * t101 / 0.2e1) * t101 + (-t78 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,2) * t84 / 0.2e1) * t84 + (-t103 * mrSges(4,1) + t94 * mrSges(4,3) + Ifges(4,4) * t101 + Ifges(4,2) * t100 / 0.2e1) * t100 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,5) * t85 + Ifges(7,6) * t84 + Ifges(7,3) * t106 / 0.2e1) * t106 + (-t119 * mrSges(5,2) + t79 * mrSges(6,2) - t81 * mrSges(5,3) - t83 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t92) * t92 + (-t119 * mrSges(5,1) + t83 * mrSges(6,1) - t80 * mrSges(6,2) - t82 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t91 + (-Ifges(5,4) + Ifges(6,5)) * t92) * t91 + (t93 * mrSges(4,1) + t81 * mrSges(5,1) - t79 * mrSges(6,1) - t94 * mrSges(4,2) - t82 * mrSges(5,2) + t80 * mrSges(6,3) + Ifges(4,5) * t101 + Ifges(4,6) * t100 + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t107 + (Ifges(6,4) + Ifges(5,5)) * t92 + (-Ifges(5,6) + Ifges(6,6)) * t91) * t107 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t112 ^ 2 + t115 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t124 + Ifges(3,2) / 0.2e1) * t115) * t115 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t115 + (t124 + Ifges(3,1) / 0.2e1) * t112) * t112) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t115 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t112) * qJD(2)) * qJD(1);
T  = t1;
