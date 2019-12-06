% Calculate potential energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:54
% EndTime: 2019-12-05 15:02:54
% DurationCPUTime: 0.38s
% Computational Cost: add. (133->65), mult. (148->62), div. (0->0), fcn. (121->8), ass. (0->30)
t87 = -mrSges(4,1) + mrSges(5,2);
t86 = mrSges(4,2) - mrSges(5,3);
t85 = m(5) + m(6);
t62 = cos(pkin(7));
t58 = pkin(8) + qJ(3);
t55 = sin(t58);
t74 = qJ(4) * t55;
t56 = cos(t58);
t78 = t62 * t56;
t84 = pkin(3) * t78 + t62 * t74;
t59 = sin(pkin(8));
t61 = cos(pkin(8));
t83 = -m(3) * pkin(1) - t61 * mrSges(3,1) + t59 * mrSges(3,2) + t86 * t55 + t87 * t56 - mrSges(2,1);
t82 = -m(3) * qJ(2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t60 = sin(pkin(7));
t81 = t60 * t56;
t64 = sin(qJ(5));
t80 = t60 * t64;
t65 = cos(qJ(5));
t79 = t60 * t65;
t77 = t62 * t64;
t76 = t62 * t65;
t54 = t61 * pkin(2) + pkin(1);
t63 = -pkin(5) - qJ(2);
t75 = t60 * t54 + t62 * t63;
t73 = t59 * pkin(2) + qJ(1);
t49 = t62 * t54;
t70 = -t60 * t63 + t49;
t69 = pkin(3) * t81 + t60 * t74 + t75;
t1 = (-m(4) * t73 - t59 * mrSges(3,1) - t61 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - t85 * (t55 * pkin(3) + t73) + (-m(2) - m(3)) * qJ(1) + (t64 * mrSges(6,1) + t65 * mrSges(6,2) + t85 * qJ(4) - t86) * t56 + (-m(6) * pkin(6) - mrSges(6,3) + t87) * t55) * g(3) + (-mrSges(1,2) - m(4) * t75 - m(5) * t69 - m(6) * (pkin(6) * t81 + t69) - (t55 * t80 - t76) * mrSges(6,1) - (t55 * t79 + t77) * mrSges(6,2) - mrSges(6,3) * t81 + (m(6) * pkin(4) - t82) * t62 + t83 * t60) * g(2) + (-mrSges(1,1) - m(4) * t70 - m(5) * (t70 + t84) - m(6) * (pkin(6) * t78 + t49 + t84) - (t55 * t77 + t79) * mrSges(6,1) - (t55 * t76 - t80) * mrSges(6,2) - mrSges(6,3) * t78 + t83 * t62 + (-m(6) * (pkin(4) - t63) + t82) * t60) * g(1);
U = t1;
