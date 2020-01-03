% Calculate Gravitation load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:48
% EndTime: 2019-12-31 20:23:51
% DurationCPUTime: 0.82s
% Computational Cost: add. (518->99), mult. (1273->146), div. (0->0), fcn. (1563->12), ass. (0->56)
t46 = sin(qJ(5));
t50 = cos(qJ(5));
t84 = m(5) + m(6);
t86 = pkin(8) * t84 - mrSges(4,2) + mrSges(5,3);
t90 = -t46 * mrSges(6,1) - t50 * mrSges(6,2) - t86;
t88 = m(6) * pkin(4) + mrSges(6,1) * t50 - mrSges(6,2) * t46 + mrSges(5,1);
t66 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t45 = cos(pkin(5));
t52 = cos(qJ(2));
t53 = cos(qJ(1));
t76 = t52 * t53;
t48 = sin(qJ(2));
t49 = sin(qJ(1));
t78 = t49 * t48;
t89 = t45 * t76 - t78;
t43 = sin(pkin(10));
t72 = cos(pkin(10));
t61 = -t43 * t48 + t52 * t72;
t47 = sin(qJ(4));
t51 = cos(qJ(4));
t87 = t47 * t66 - t51 * t88 - mrSges(4,1);
t83 = pkin(2) * t52;
t44 = sin(pkin(5));
t82 = t44 * t49;
t81 = t44 * t53;
t79 = t48 * t53;
t77 = t49 * t52;
t33 = -t43 * t52 - t48 * t72;
t26 = t33 * t45;
t16 = -t26 * t49 - t53 * t61;
t42 = pkin(1) + t83;
t36 = t53 * t42;
t74 = -pkin(3) * t16 + t36;
t71 = m(4) + t84;
t40 = t44 * t83;
t69 = -m(3) * pkin(1) - mrSges(2,1);
t15 = -t26 * t53 + t49 * t61;
t4 = t15 * t51 - t47 * t81;
t3 = -t15 * t47 - t51 * t81;
t65 = t89 * pkin(2);
t30 = -t45 * t77 - t79;
t58 = t30 * pkin(2);
t57 = t45 * t61;
t55 = mrSges(2,2) + (-m(3) * pkin(7) - mrSges(3,3) - mrSges(4,3)) * t44 + t71 * (pkin(2) * t45 * t48 + (-pkin(7) - qJ(3)) * t44);
t31 = -t45 * t78 + t76;
t29 = -t45 * t79 - t77;
t25 = t33 * t44;
t24 = t61 * t44;
t20 = -t25 * t51 + t45 * t47;
t17 = t33 * t53 - t49 * t57;
t14 = t49 * t33 + t53 * t57;
t8 = -t16 * t51 + t47 * t82;
t7 = -t16 * t47 - t51 * t82;
t2 = -t17 * t46 + t50 * t8;
t1 = -t17 * t50 - t46 * t8;
t5 = [(-t31 * mrSges(3,1) - t30 * mrSges(3,2) - m(4) * t36 + t16 * mrSges(4,1) - m(5) * t74 - t8 * mrSges(5,1) - m(6) * (pkin(4) * t8 + t74) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t66 * t7 + t69 * t53 + t86 * t17 + t55 * t49) * g(2) + (-t29 * mrSges(3,1) + t89 * mrSges(3,2) + t66 * t3 + t88 * t4 + t90 * t14 + (t42 * t71 - t69) * t49 + t55 * t53 + (pkin(3) * t84 + mrSges(4,1)) * t15) * g(1), (-(mrSges(3,1) * t52 - mrSges(3,2) * t48) * t44 - m(4) * t40 - t84 * (pkin(3) * t24 + t40) - t90 * t25 + t87 * t24) * g(3) + (-m(4) * t65 - t89 * mrSges(3,1) - t29 * mrSges(3,2) - t84 * (pkin(3) * t14 + t65) + t87 * t14 + t90 * t15) * g(2) + (-m(4) * t58 - t30 * mrSges(3,1) + t31 * mrSges(3,2) - t84 * (pkin(3) * t17 + t58) + t87 * t17 - t90 * t16) * g(1), (-g(3) * t45 + (-g(1) * t49 + g(2) * t53) * t44) * t71, (t66 * t20 - t88 * (t25 * t47 + t45 * t51)) * g(3) + (-t3 * t88 + t4 * t66) * g(2) + (t66 * t8 + t7 * t88) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((-t14 * t50 - t4 * t46) * mrSges(6,1) + (t14 * t46 - t4 * t50) * mrSges(6,2)) - g(3) * ((-t20 * t46 - t24 * t50) * mrSges(6,1) + (-t20 * t50 + t24 * t46) * mrSges(6,2))];
taug = t5(:);
