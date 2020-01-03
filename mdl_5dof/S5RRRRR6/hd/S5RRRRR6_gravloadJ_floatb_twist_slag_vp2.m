% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:47
% EndTime: 2020-01-03 12:14:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (327->77), mult. (253->83), div. (0->0), fcn. (194->10), ass. (0->51)
t43 = qJ(3) + qJ(4);
t38 = qJ(5) + t43;
t32 = cos(t38);
t19 = t32 * mrSges(6,1);
t36 = cos(t43);
t76 = -t36 * mrSges(5,1) - t19;
t47 = cos(qJ(3));
t45 = sin(qJ(3));
t34 = sin(t43);
t31 = sin(t38);
t68 = t31 * mrSges(6,2);
t54 = t34 * mrSges(5,2) + t68;
t52 = t45 * mrSges(4,2) + t54;
t78 = -t47 * mrSges(4,1) - mrSges(3,1) + t52;
t77 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t60 = m(5) * pkin(3) + mrSges(4,1);
t75 = -t60 * t45 + m(6) * (-t45 * pkin(3) - pkin(4) * t34) - mrSges(4,2) * t47;
t74 = m(6) * pkin(4);
t49 = -pkin(8) - pkin(7);
t44 = qJ(1) + qJ(2);
t37 = cos(t44);
t67 = t32 * t37;
t69 = mrSges(6,1) * t31;
t73 = mrSges(6,2) * t67 + t37 * t69;
t30 = pkin(4) * t36;
t35 = sin(t44);
t72 = g(2) * t35;
t42 = -pkin(9) + t49;
t40 = t47 * pkin(3);
t62 = t30 + t40;
t7 = pkin(2) + t62;
t71 = t35 * t7 + t37 * t42;
t33 = t40 + pkin(2);
t70 = t35 * t33 + t37 * t49;
t66 = t36 * t37;
t65 = t37 * t34;
t63 = t37 * pkin(2) + t35 * pkin(7);
t61 = -m(3) * pkin(1) - mrSges(2,1);
t59 = -t35 * t42 + t37 * t7;
t58 = -mrSges(5,1) * t65 - mrSges(5,2) * t66 - t73;
t57 = t37 * t33 - t35 * t49;
t55 = -mrSges(6,2) * t32 - t69;
t53 = mrSges(5,2) * t36 - t55;
t51 = -mrSges(5,1) * t66 - mrSges(6,1) * t67 + t77 * t35 + t78 * t37;
t50 = (m(4) * pkin(7) - t77) * t37 + (t78 + t76) * t35;
t48 = cos(qJ(1));
t46 = sin(qJ(1));
t41 = t48 * pkin(1);
t39 = t46 * pkin(1);
t28 = t35 * pkin(2);
t1 = [(-t48 * mrSges(2,2) - m(4) * (t28 + t39) - m(5) * (t39 + t70) - m(6) * (t39 + t71) + t61 * t46 + t50) * g(3) + (t46 * mrSges(2,2) - m(4) * (t41 + t63) - m(5) * (t41 + t57) - m(6) * (t41 + t59) + t61 * t48 + t51) * g(2), (-m(4) * t28 - m(5) * t70 - m(6) * t71 + t50) * g(3) + (-m(4) * t63 - m(5) * t57 - m(6) * t59 + t51) * g(2), (t75 * t37 + t58) * g(3) + (-m(6) * t62 - t60 * t47 + t52 + t76) * g(1) + (mrSges(5,1) * t34 + t53 - t75) * t72, (-t65 * t74 + t58) * g(3) + (-m(6) * t30 + t54 + t76) * g(1) + ((mrSges(5,1) + t74) * t34 + t53) * t72, -g(1) * (t19 - t68) - g(3) * t73 - t55 * t72];
taug = t1(:);
