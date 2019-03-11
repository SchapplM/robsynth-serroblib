% Calculate kinetic energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:41
% EndTime: 2019-03-09 02:46:42
% DurationCPUTime: 0.56s
% Computational Cost: add. (547->101), mult. (1337->145), div. (0->0), fcn. (958->8), ass. (0->43)
t113 = cos(pkin(9));
t130 = t113 ^ 2;
t129 = m(3) / 0.2e1;
t128 = -pkin(4) - pkin(5);
t127 = cos(qJ(3));
t126 = pkin(7) + qJ(2);
t111 = sin(pkin(10));
t125 = cos(pkin(10));
t115 = sin(qJ(3));
t122 = qJD(1) * t113;
t112 = sin(pkin(9));
t123 = qJD(1) * t112;
t101 = t115 * t123 - t127 * t122;
t102 = (t127 * t112 + t113 * t115) * qJD(1);
t105 = qJD(2) + (-pkin(2) * t113 - pkin(1)) * qJD(1);
t86 = pkin(3) * t101 - qJ(4) * t102 + t105;
t103 = t126 * t123;
t104 = t126 * t122;
t92 = -t115 * t103 + t127 * t104;
t90 = qJD(3) * qJ(4) + t92;
t82 = t111 * t86 + t125 * t90;
t91 = -t127 * t103 - t115 * t104;
t79 = t101 * qJ(5) + t82;
t81 = -t111 * t90 + t125 * t86;
t121 = qJD(3) * pkin(3) - qJD(4) + t91;
t120 = qJD(5) - t81;
t94 = t111 * qJD(3) + t125 * t102;
t119 = qJ(5) * t94 + t121;
t116 = cos(qJ(6));
t114 = sin(qJ(6));
t107 = -qJD(1) * pkin(1) + qJD(2);
t96 = qJD(6) - t101;
t93 = -t125 * qJD(3) + t102 * t111;
t84 = t114 * t93 + t116 * t94;
t83 = -t114 * t94 + t116 * t93;
t80 = pkin(4) * t93 - t119;
t78 = -t101 * pkin(4) + t120;
t77 = t128 * t93 + t119;
t76 = pkin(8) * t93 + t79;
t75 = -t94 * pkin(8) + t128 * t101 + t120;
t74 = t114 * t75 + t116 * t76;
t73 = -t114 * t76 + t116 * t75;
t1 = t107 ^ 2 * t129 + m(4) * (t105 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + (t73 * mrSges(7,1) - t74 * mrSges(7,2) + Ifges(7,3) * t96 / 0.2e1) * t96 + (t105 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,1) * t102 / 0.2e1) * t102 + (t77 * mrSges(7,2) - t73 * mrSges(7,3) + Ifges(7,5) * t96 + Ifges(7,1) * t84 / 0.2e1) * t84 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,5) * t102 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t77 * mrSges(7,1) + t74 * mrSges(7,3) + Ifges(7,4) * t84 + Ifges(7,6) * t96 + Ifges(7,2) * t83 / 0.2e1) * t83 + (-t121 * mrSges(5,2) + t78 * mrSges(6,2) - t81 * mrSges(5,3) - t80 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t94) * t94 + (-t121 * mrSges(5,1) + t80 * mrSges(6,1) - t79 * mrSges(6,2) - t82 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t93 + (-Ifges(5,4) + Ifges(6,5)) * t94) * t93 + (t105 * mrSges(4,1) + t81 * mrSges(5,1) - t78 * mrSges(6,1) - t82 * mrSges(5,2) - t92 * mrSges(4,3) + t79 * mrSges(6,3) - Ifges(4,4) * t102 - Ifges(4,6) * qJD(3) + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101 + (Ifges(6,4) + Ifges(5,5)) * t94 + (-Ifges(5,6) + Ifges(6,6)) * t93) * t101 + (t107 * (-mrSges(3,1) * t113 + mrSges(3,2) * t112) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t129 + mrSges(3,3)) * (t112 ^ 2 + t130) * qJ(2) + Ifges(3,2) * t130 / 0.2e1 + (Ifges(3,4) * t113 + Ifges(3,1) * t112 / 0.2e1) * t112) * qJD(1)) * qJD(1);
T  = t1;
