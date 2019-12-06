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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 19:00:03
% EndTime: 2019-12-05 19:00:03
% DurationCPUTime: 0.23s
% Computational Cost: add. (327->64), mult. (253->62), div. (0->0), fcn. (194->10), ass. (0->42)
t33 = qJ(3) + qJ(4);
t30 = qJ(5) + t33;
t24 = cos(t30);
t15 = t24 * mrSges(6,1);
t28 = cos(t33);
t17 = t28 * mrSges(5,1);
t37 = cos(qJ(3));
t31 = t37 * pkin(3);
t22 = pkin(4) * t28;
t51 = t22 + t31;
t35 = sin(qJ(3));
t52 = t35 * mrSges(4,2);
t67 = -t52 + m(4) * pkin(2) + m(5) * (t31 + pkin(2)) + m(6) * (pkin(2) + t51) + t37 * mrSges(4,1) + mrSges(3,1) + t15 + t17;
t39 = -pkin(8) - pkin(7);
t66 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-pkin(9) + t39) + m(5) * t39;
t26 = sin(t33);
t49 = m(5) * pkin(3) + mrSges(4,1);
t62 = -t49 * t35 + m(6) * (-t35 * pkin(3) - pkin(4) * t26) - mrSges(4,2) * t37;
t61 = m(6) * pkin(4);
t34 = qJ(1) + qJ(2);
t27 = sin(t34);
t23 = sin(t30);
t55 = t23 * t27;
t57 = mrSges(6,2) * t24;
t60 = mrSges(6,1) * t55 + t27 * t57;
t29 = cos(t34);
t59 = g(3) * t29;
t58 = mrSges(5,2) * t28;
t56 = t23 * mrSges(6,2);
t54 = t26 * mrSges(5,2);
t53 = t26 * t27;
t50 = -mrSges(5,1) * t53 - t27 * t58 - t60;
t48 = t15 - t56;
t46 = -mrSges(6,1) * t23 - t57;
t45 = mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1);
t44 = -t17 - t48 + t54;
t43 = -t46 + t58;
t41 = (-t54 - t56 + t67) * t29 - t66 * t27;
t40 = -mrSges(5,2) * t53 - mrSges(6,2) * t55 + t67 * t27 + t66 * t29;
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t1 = [(t38 * mrSges(2,2) + t45 * t36 + t40) * g(3) + (-t36 * mrSges(2,2) + t45 * t38 + t41) * g(2), t41 * g(2) + t40 * g(3), (t62 * t27 + t50) * g(2) + (-m(6) * t51 - t49 * t37 + t44 + t52) * g(1) + (mrSges(5,1) * t26 + t43 - t62) * t59, (-t53 * t61 + t50) * g(2) + (-m(6) * t22 + t44) * g(1) + ((mrSges(5,1) + t61) * t26 + t43) * t59, -g(1) * t48 - g(2) * t60 - t46 * t59];
taug = t1(:);
