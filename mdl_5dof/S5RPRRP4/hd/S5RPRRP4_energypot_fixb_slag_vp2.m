% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:09
% EndTime: 2020-01-03 11:49:09
% DurationCPUTime: 0.51s
% Computational Cost: add. (142->70), mult. (187->65), div. (0->0), fcn. (168->8), ass. (0->29)
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t82 = -m(4) * pkin(2) - t61 * mrSges(4,1) + t59 * mrSges(4,2) - mrSges(3,1);
t81 = mrSges(3,2) - mrSges(5,3) - mrSges(6,3);
t80 = -mrSges(5,1) - mrSges(6,1);
t79 = mrSges(5,2) + mrSges(6,2);
t78 = -m(3) - m(5) - m(6);
t77 = m(4) * qJ(2) + mrSges(4,1) * t59 + mrSges(4,2) * t61 - mrSges(2,2) + mrSges(3,3);
t76 = m(4) * pkin(6) + mrSges(4,3);
t57 = sin(pkin(8));
t58 = cos(pkin(8));
t75 = t81 * t57 + t82 * t58 - mrSges(2,1);
t63 = -pkin(7) - pkin(6);
t74 = pkin(3) * t59;
t50 = t61 * pkin(3) + pkin(2);
t60 = sin(qJ(1));
t71 = t58 * t60;
t62 = cos(qJ(1));
t70 = t58 * t62;
t56 = qJ(3) + qJ(4);
t52 = cos(t56);
t48 = pkin(4) * t52 + t50;
t55 = -qJ(5) + t63;
t67 = t48 * t58 - t55 * t57;
t66 = t50 * t58 - t57 * t63;
t53 = t60 * pkin(1);
t51 = sin(t56);
t49 = pkin(4) * t51 + t74;
t1 = (-mrSges(1,3) + t80 * (-t51 * t60 - t52 * t70) - t79 * (t51 * t70 - t52 * t60) + (m(3) * pkin(1) - m(4) * (-pkin(6) * t57 - pkin(1)) + t57 * mrSges(4,3) - m(5) * (-pkin(1) - t66) - m(6) * (-pkin(1) - t67) - t75) * t62 + (m(3) * qJ(2) - m(5) * (-qJ(2) - t74) - m(6) * (-qJ(2) - t49) + t77) * t60) * g(3) + (-m(4) * t53 - mrSges(1,2) + t80 * (-t51 * t62 + t52 * t71) - t79 * (-t51 * t71 - t52 * t62) + t78 * (-t62 * qJ(2) + t53) + (m(5) * t74 + m(6) * t49 + t77) * t62 + (-m(5) * t66 - m(6) * t67 - t76 * t57 + t75) * t60) * g(2) + (-mrSges(1,1) - mrSges(2,3) + (-m(2) - m(4) + t78) * pkin(5) + (-m(5) * t63 - m(6) * t55 + t76 - t81) * t58 + (-m(5) * t50 - m(6) * t48 + t79 * t51 + t80 * t52 + t82) * t57) * g(1);
U = t1;
