% Calculate potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:31
% EndTime: 2019-12-05 17:35:32
% DurationCPUTime: 0.38s
% Computational Cost: add. (147->60), mult. (148->59), div. (0->0), fcn. (125->8), ass. (0->27)
t81 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t80 = m(6) * pkin(4);
t79 = -mrSges(5,1) - mrSges(6,1);
t78 = mrSges(3,2) - mrSges(4,3);
t77 = mrSges(5,2) + mrSges(6,2);
t76 = -m(4) - m(5) - m(6);
t55 = sin(pkin(8));
t56 = cos(pkin(8));
t75 = t56 * mrSges(4,1) + t81 * t55 + mrSges(3,1);
t60 = sin(qJ(1));
t74 = pkin(1) * t60;
t62 = cos(qJ(1));
t53 = t62 * pkin(1);
t54 = qJ(1) + pkin(7);
t51 = sin(t54);
t59 = sin(qJ(4));
t73 = t51 * t59;
t52 = cos(t54);
t72 = t52 * t59;
t69 = t56 * t59;
t61 = cos(qJ(4));
t68 = t56 * t61;
t65 = pkin(3) * t56 + pkin(6) * t55;
t50 = pkin(4) * t61 + pkin(3);
t57 = -qJ(5) - pkin(6);
t63 = t50 * t56 - t55 * t57;
t1 = (-t73 * t80 - m(3) * t53 - mrSges(2,1) * t62 + t60 * mrSges(2,2) - mrSges(1,3) + t76 * (t52 * pkin(2) + t51 * qJ(3) + t53) + t78 * t51 + t79 * (t52 * t68 + t73) - t77 * (t51 * t61 - t52 * t69) + (-m(5) * t65 - m(6) * t63 - t75) * t52) * g(3) + (-t72 * t80 + m(3) * t74 + t60 * mrSges(2,1) + mrSges(2,2) * t62 - mrSges(1,2) + t76 * (t52 * qJ(3) - t74) + t78 * t52 + t79 * (-t51 * t68 + t72) - t77 * (t51 * t69 + t52 * t61) + (m(4) * pkin(2) - m(5) * (-pkin(2) - t65) - m(6) * (-pkin(2) - t63) + t75) * t51) * g(2) + (-m(2) * pkin(5) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t76) * (qJ(2) + pkin(5)) + (m(5) * pkin(6) - m(6) * t57 + t81) * t56 + (-m(5) * pkin(3) - m(6) * t50 + t77 * t59 + t79 * t61 - mrSges(4,1)) * t55) * g(1);
U = t1;
