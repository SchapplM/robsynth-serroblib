% Calculate Gravitation load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:02
% EndTime: 2018-11-23 15:15:03
% DurationCPUTime: 1.13s
% Computational Cost: add. (1204->107), mult. (1229->146), div. (0->0), fcn. (1182->18), ass. (0->50)
t43 = qJ(5) + qJ(6);
t40 = sin(t43);
t41 = cos(t43);
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t114 = mrSges(5,1) + m(7) * (pkin(5) * t51 + pkin(4)) + t41 * mrSges(7,1) - t40 * mrSges(7,2) + m(6) * pkin(4) + t51 * mrSges(6,1) - t48 * mrSges(6,2);
t99 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t42 = qJ(3) + pkin(12);
t38 = sin(t42);
t39 = cos(t42);
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t101 = m(4) * pkin(2) + t52 * mrSges(4,1) - t49 * mrSges(4,2) + t114 * t39 - t99 * t38 + mrSges(3,1);
t108 = -m(7) * pkin(5) - mrSges(6,1);
t104 = m(5) + m(6) + m(7);
t84 = pkin(6) + qJ(2);
t71 = sin(t84) / 0.2e1;
t85 = pkin(6) - qJ(2);
t77 = sin(t85);
t28 = t71 - t77 / 0.2e1;
t44 = sin(pkin(11));
t53 = cos(qJ(2));
t86 = cos(pkin(11));
t65 = -t28 * t44 + t53 * t86;
t45 = sin(pkin(6));
t91 = t44 * t45;
t111 = -t49 * t65 + t52 * t91;
t72 = cos(t84) / 0.2e1;
t78 = cos(t85);
t29 = t72 - t78 / 0.2e1;
t46 = cos(pkin(6));
t110 = t29 * t49 + t46 * t52;
t66 = t28 * t86 + t44 * t53;
t79 = t45 * t86;
t106 = -t49 * t66 - t52 * t79;
t100 = -m(4) * pkin(8) - t40 * mrSges(7,1) - t51 * mrSges(6,2) - t41 * mrSges(7,2) + t108 * t48 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t50 = sin(qJ(2));
t57 = t78 / 0.2e1 + t72;
t17 = t44 * t50 - t57 * t86;
t8 = -t38 * t79 + t39 * t66;
t97 = (t17 * t41 - t40 * t8) * mrSges(7,1) + (-t17 * t40 - t41 * t8) * mrSges(7,2);
t10 = t38 * t91 + t39 * t65;
t20 = t44 * t57 + t50 * t86;
t96 = (-t10 * t40 + t20 * t41) * mrSges(7,1) + (-t10 * t41 - t20 * t40) * mrSges(7,2);
t14 = -t29 * t39 + t38 * t46;
t27 = t71 + t77 / 0.2e1;
t95 = (-t14 * t40 - t27 * t41) * mrSges(7,1) + (-t14 * t41 + t27 * t40) * mrSges(7,2);
t47 = -qJ(4) - pkin(8);
t37 = pkin(3) * t52 + pkin(2);
t1 = [(-m(2) - m(3) - m(4) - t104) * g(3) (-t104 * (t27 * t37 + t29 * t47) - t100 * t29 - t101 * t27) * g(3) + (-t104 * (-t17 * t37 - t47 * t66) + t100 * t66 + t101 * t17) * g(2) + (-t104 * (-t20 * t37 - t47 * t65) + t100 * t65 + t101 * t20) * g(1) (-t110 * mrSges(4,1) - (t29 * t52 - t46 * t49) * mrSges(4,2) + t99 * t14 - t114 * (t29 * t38 + t39 * t46)) * g(3) + (-(t49 * t79 - t52 * t66) * mrSges(4,2) - mrSges(4,1) * t106 + t99 * t8 - t114 * (-t38 * t66 - t39 * t79)) * g(2) + (-t111 * mrSges(4,1) - (-t49 * t91 - t52 * t65) * mrSges(4,2) - t114 * (-t38 * t65 + t39 * t91) + t99 * t10) * g(1) + (-g(1) * t111 - t106 * g(2) - g(3) * t110) * t104 * pkin(3), t104 * (-g(1) * t20 - g(2) * t17 + g(3) * t27) (-(-t14 * t51 + t27 * t48) * mrSges(6,2) - t95 + t108 * (-t14 * t48 - t27 * t51)) * g(3) + (-(-t17 * t48 - t51 * t8) * mrSges(6,2) - t97 + t108 * (t17 * t51 - t48 * t8)) * g(2) + (-(-t10 * t51 - t20 * t48) * mrSges(6,2) - t96 + t108 * (-t10 * t48 + t20 * t51)) * g(1), -g(1) * t96 - g(2) * t97 - g(3) * t95];
taug  = t1(:);
