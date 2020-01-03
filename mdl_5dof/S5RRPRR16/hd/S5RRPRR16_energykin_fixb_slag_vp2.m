% Calculate kinetic energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR16_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:38
% EndTime: 2019-12-31 20:44:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (376->89), mult. (867->136), div. (0->0), fcn. (590->8), ass. (0->40)
t111 = -pkin(2) - pkin(8);
t102 = cos(qJ(4));
t100 = sin(qJ(2));
t103 = cos(qJ(2));
t96 = sin(pkin(5));
t110 = qJD(1) * t96;
t107 = t100 * t110;
t91 = pkin(7) * t107;
t97 = cos(pkin(5));
t109 = qJD(1) * t97;
t95 = qJD(2) + t109;
t75 = qJD(3) + t91 + t111 * t95 + (-pkin(1) * t103 * t97 + pkin(3) * t100 * t96) * qJD(1);
t105 = -qJ(3) * t100 - pkin(1);
t80 = (t111 * t103 + t105) * t110;
t99 = sin(qJ(4));
t72 = t102 * t80 + t99 * t75;
t106 = t103 * t110;
t108 = pkin(1) * t109;
t88 = pkin(7) * t106 + t100 * t108;
t82 = -t95 * qJ(3) - t88;
t79 = pkin(3) * t106 - t82;
t71 = t102 * t75 - t80 * t99;
t87 = t103 * t108 - t91;
t85 = -t102 * t106 - t95 * t99;
t104 = qJD(1) ^ 2;
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t89 = qJD(4) + t107;
t86 = t102 * t95 - t99 * t106;
t84 = qJD(5) - t85;
t83 = (-pkin(2) * t103 + t105) * t110;
t81 = -pkin(2) * t95 + qJD(3) - t87;
t77 = t101 * t86 + t89 * t98;
t76 = t101 * t89 - t86 * t98;
t73 = -pkin(4) * t85 - pkin(9) * t86 + t79;
t70 = pkin(9) * t89 + t72;
t69 = -pkin(4) * t89 - t71;
t68 = t101 * t70 + t73 * t98;
t67 = t101 * t73 - t70 * t98;
t1 = m(3) * (pkin(1) ^ 2 * t104 * t96 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t104 * Ifges(2,3) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + (t71 * mrSges(5,1) - t72 * mrSges(5,2) + Ifges(5,3) * t89 / 0.2e1) * t89 + (t67 * mrSges(6,1) - t68 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t79 * mrSges(5,2) - t71 * mrSges(5,3) + Ifges(5,5) * t89 + Ifges(5,1) * t86 / 0.2e1) * t86 + (t69 * mrSges(6,2) - t67 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t77 / 0.2e1) * t77 + (-t79 * mrSges(5,1) + t72 * mrSges(5,3) + Ifges(5,4) * t86 + Ifges(5,6) * t89 + Ifges(5,2) * t85 / 0.2e1) * t85 + (-t69 * mrSges(6,1) + t68 * mrSges(6,3) + Ifges(6,4) * t77 + Ifges(6,6) * t84 + Ifges(6,2) * t76 / 0.2e1) * t76 + ((-t82 * mrSges(4,1) + t83 * mrSges(4,2) + t88 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t103) * t110) * t103 + (t81 * mrSges(4,1) - t87 * mrSges(3,3) - t83 * mrSges(4,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t100 + (Ifges(3,4) + Ifges(4,6)) * t103) * t110) * t100) * t110 + (t87 * mrSges(3,1) - t88 * mrSges(3,2) + t81 * mrSges(4,2) - t82 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t95 + ((-Ifges(4,5) + Ifges(3,6)) * t103 + (-Ifges(4,4) + Ifges(3,5)) * t100) * t110) * t95;
T = t1;
