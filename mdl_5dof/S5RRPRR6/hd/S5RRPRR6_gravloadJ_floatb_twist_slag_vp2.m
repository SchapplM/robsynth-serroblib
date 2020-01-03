% Calculate Gravitation load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:29
% EndTime: 2020-01-03 12:05:31
% DurationCPUTime: 0.33s
% Computational Cost: add. (334->83), mult. (316->101), div. (0->0), fcn. (294->10), ass. (0->52)
t84 = mrSges(3,2) - mrSges(4,3);
t83 = -mrSges(6,3) - mrSges(5,3);
t50 = sin(pkin(9));
t82 = t50 * mrSges(4,2) - mrSges(3,1);
t80 = m(6) * pkin(4);
t81 = -mrSges(5,1) - t80;
t48 = qJ(4) + qJ(5);
t42 = sin(t48);
t44 = cos(t48);
t49 = qJ(1) + qJ(2);
t45 = cos(t49);
t43 = sin(t49);
t51 = cos(pkin(9));
t75 = t43 * t51;
t10 = -t45 * t42 + t44 * t75;
t9 = -t42 * t75 - t45 * t44;
t79 = t9 * mrSges(6,1) - t10 * mrSges(6,2);
t72 = t45 * t51;
t11 = t42 * t72 - t43 * t44;
t12 = t43 * t42 + t44 * t72;
t78 = t11 * mrSges(6,1) + t12 * mrSges(6,2);
t77 = g(1) * t50;
t76 = t43 * t50;
t52 = sin(qJ(4));
t74 = t43 * t52;
t73 = t45 * t50;
t70 = t50 * (-pkin(8) - pkin(7));
t69 = t51 * t52;
t54 = cos(qJ(4));
t68 = t51 * t54;
t67 = t45 * pkin(2) + t43 * qJ(3);
t66 = m(4) + m(5) + m(6);
t65 = t52 * t80;
t39 = t43 * pkin(2);
t64 = pkin(3) * t75 + pkin(7) * t76 + t39;
t63 = -m(3) * pkin(1) - mrSges(2,1);
t62 = pkin(3) * t72 + pkin(7) * t73 + t67;
t61 = -mrSges(6,1) * t42 - mrSges(6,2) * t44;
t41 = t54 * pkin(4) + pkin(3);
t60 = t41 * t75 - t43 * t70 + t39;
t19 = -t43 * t54 + t45 * t69;
t17 = -t43 * t69 - t45 * t54;
t59 = pkin(4) * t74 + t41 * t72 - t45 * t70 + t67;
t20 = t45 * t68 + t74;
t58 = -mrSges(4,1) * t72 - t20 * mrSges(5,1) - t12 * mrSges(6,1) + t19 * mrSges(5,2) + t11 * mrSges(6,2) + t84 * t43 + t82 * t45 + t83 * t73;
t18 = t43 * t68 - t45 * t52;
t57 = -mrSges(4,1) * t75 - t18 * mrSges(5,1) - t10 * mrSges(6,1) - t17 * mrSges(5,2) - t9 * mrSges(6,2) + t83 * t76 + (t66 * qJ(3) + t65 - t84) * t45 + t82 * t43;
t55 = cos(qJ(1));
t53 = sin(qJ(1));
t47 = t55 * pkin(1);
t46 = t53 * pkin(1);
t1 = [(t57 + t63 * t53 - m(5) * (t46 + t64) - m(4) * (t39 + t46) - m(6) * (t46 + t60) - t55 * mrSges(2,2)) * g(3) + (t58 + t63 * t55 - m(4) * (t47 + t67) - m(5) * (t47 + t62) + t53 * mrSges(2,2) - m(6) * (t47 + t59)) * g(2), (-m(4) * t39 - m(5) * t64 - m(6) * t60 + t57) * g(3) + (-m(4) * t67 - m(5) * t62 - m(6) * t59 + t58) * g(2), (g(2) * t45 + g(3) * t43) * t66, (mrSges(5,1) * t52 + mrSges(5,2) * t54 - t61 + t65) * t77 + (-t20 * mrSges(5,2) + t81 * t19 - t78) * g(3) + (t18 * mrSges(5,2) + t81 * t17 - t79) * g(2), -g(2) * t79 - g(3) * t78 - t61 * t77];
taug = t1(:);
