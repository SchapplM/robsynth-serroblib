% Calculate kinetic energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:39
% EndTime: 2019-12-31 19:09:40
% DurationCPUTime: 0.45s
% Computational Cost: add. (403->80), mult. (970->129), div. (0->0), fcn. (686->8), ass. (0->38)
t98 = cos(pkin(9));
t112 = t98 ^ 2;
t111 = m(3) / 0.2e1;
t110 = pkin(6) + qJ(2);
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t101 = sin(qJ(3));
t104 = cos(qJ(3));
t107 = qJD(1) * t98;
t97 = sin(pkin(9));
t108 = qJD(1) * t97;
t88 = -t101 * t108 + t104 * t107;
t89 = (t101 * t98 + t104 * t97) * qJD(1);
t92 = qJD(2) + (-pkin(2) * t98 - pkin(1)) * qJD(1);
t75 = -pkin(3) * t88 - pkin(7) * t89 + t92;
t90 = t110 * t108;
t91 = t110 * t107;
t80 = -t101 * t90 + t104 * t91;
t78 = qJD(3) * pkin(7) + t80;
t69 = t100 * t75 + t103 * t78;
t68 = -t100 * t78 + t103 * t75;
t79 = -t101 * t91 - t104 * t90;
t84 = qJD(4) - t88;
t77 = -qJD(3) * pkin(3) - t79;
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t94 = -qJD(1) * pkin(1) + qJD(2);
t83 = qJD(5) + t84;
t82 = qJD(3) * t100 + t103 * t89;
t81 = qJD(3) * t103 - t100 * t89;
t72 = t102 * t82 + t81 * t99;
t71 = t102 * t81 - t82 * t99;
t70 = -pkin(4) * t81 + t77;
t67 = pkin(8) * t81 + t69;
t66 = pkin(4) * t84 - pkin(8) * t82 + t68;
t65 = t102 * t67 + t66 * t99;
t64 = t102 * t66 - t67 * t99;
t1 = m(4) * (t79 ^ 2 + t80 ^ 2 + t92 ^ 2) / 0.2e1 + t94 ^ 2 * t111 + m(6) * (t64 ^ 2 + t65 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t77 ^ 2) / 0.2e1 + (t92 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,1) * t89 / 0.2e1) * t89 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + Ifges(5,3) * t84 / 0.2e1) * t84 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (-t92 * mrSges(4,1) + t80 * mrSges(4,3) + Ifges(4,4) * t89 + Ifges(4,2) * t88 / 0.2e1) * t88 + (t77 * mrSges(5,2) - t68 * mrSges(5,3) + Ifges(5,5) * t84 + Ifges(5,1) * t82 / 0.2e1) * t82 + (t70 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t72 / 0.2e1) * t72 + (-t77 * mrSges(5,1) + t69 * mrSges(5,3) + Ifges(5,4) * t82 + Ifges(5,6) * t84 + Ifges(5,2) * t81 / 0.2e1) * t81 + (-t70 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t72 + Ifges(6,6) * t83 + Ifges(6,2) * t71 / 0.2e1) * t71 + (t79 * mrSges(4,1) - t80 * mrSges(4,2) + Ifges(4,5) * t89 + Ifges(4,6) * t88 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t94 * (-mrSges(3,1) * t98 + mrSges(3,2) * t97) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t111 + mrSges(3,3)) * (t97 ^ 2 + t112) * qJ(2) + Ifges(3,2) * t112 / 0.2e1 + (Ifges(3,4) * t98 + Ifges(3,1) * t97 / 0.2e1) * t97) * qJD(1)) * qJD(1);
T = t1;
