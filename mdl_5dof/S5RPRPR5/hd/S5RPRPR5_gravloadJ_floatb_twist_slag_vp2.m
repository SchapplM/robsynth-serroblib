% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:56
% EndTime: 2022-01-23 09:24:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (233->93), mult. (299->112), div. (0->0), fcn. (277->10), ass. (0->57)
t40 = cos(qJ(3));
t32 = t40 * pkin(3);
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t48 = qJ(4) + pkin(6);
t66 = pkin(2) * t37 + pkin(1);
t79 = m(4) * t66 + m(5) * (t37 * t32 + t66) + mrSges(2,1) + t37 * mrSges(3,1) + (m(4) * pkin(6) + m(5) * t48 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t36;
t77 = -m(3) - m(4);
t71 = m(5) + m(6);
t76 = m(5) * pkin(3) + mrSges(4,1);
t38 = sin(qJ(3));
t67 = t38 * pkin(3);
t74 = -m(5) * (qJ(2) + t67) + mrSges(2,2) - mrSges(3,3);
t35 = qJ(3) + pkin(9);
t28 = qJ(5) + t35;
t24 = cos(t28);
t41 = cos(qJ(1));
t54 = t41 * t24;
t23 = sin(t28);
t39 = sin(qJ(1));
t63 = t39 * t23;
t5 = t37 * t63 + t54;
t55 = t41 * t23;
t62 = t39 * t24;
t6 = -t37 * t62 + t55;
t70 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t37 * t55 + t62;
t8 = t37 * t54 + t63;
t69 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t68 = g(3) * t36;
t26 = sin(t35);
t19 = pkin(4) * t26 + t67;
t64 = t39 * t19;
t61 = t39 * t26;
t27 = cos(t35);
t60 = t39 * t27;
t58 = t39 * t38;
t57 = t39 * t40;
t56 = t41 * t19;
t53 = t41 * t26;
t52 = t41 * t27;
t50 = t41 * t38;
t49 = t41 * t40;
t20 = pkin(4) * t27 + t32;
t29 = t39 * qJ(2);
t47 = t41 * pkin(1) + t29;
t43 = -mrSges(6,1) * t23 - mrSges(6,2) * t24;
t42 = (pkin(2) + t20) * t37 - (-pkin(7) - t48) * t36;
t16 = -t37 * t50 + t57;
t14 = t37 * t58 + t49;
t17 = t37 * t49 + t58;
t15 = -t37 * t57 + t50;
t12 = t37 * t52 + t61;
t11 = -t37 * t53 + t60;
t10 = -t37 * t60 + t53;
t9 = t37 * t61 + t52;
t1 = [(-m(3) * t47 - m(4) * t29 - t17 * mrSges(4,1) - t16 * mrSges(4,2) - t12 * mrSges(5,1) - t11 * mrSges(5,2) - m(6) * (t47 + t64) - t8 * mrSges(6,1) - t7 * mrSges(6,2) + t74 * t39) * g(2) + (-m(6) * t56 - t15 * mrSges(4,1) - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t14 * mrSges(4,2) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(6) * (-pkin(1) - t42) + t79) * t39) * g(1) + ((-m(6) * t42 - t79) * g(2) + ((-m(6) + t77) * qJ(2) + t74) * g(1)) * t41, (-g(1) * t39 + g(2) * t41) * (t71 - t77), (m(5) * t67 + m(6) * t19 + mrSges(4,1) * t38 + mrSges(5,1) * t26 + mrSges(4,2) * t40 + mrSges(5,2) * t27 - t43) * t68 + (-t15 * mrSges(4,2) + t9 * mrSges(5,1) - t10 * mrSges(5,2) - m(6) * (-t41 * t20 - t37 * t64) - t70 + t76 * t14) * g(2) + (t17 * mrSges(4,2) - t11 * mrSges(5,1) + t12 * mrSges(5,2) - m(6) * (t39 * t20 - t37 * t56) - t69 - t76 * t16) * g(1), (t37 * g(3) + t36 * (-g(1) * t41 - g(2) * t39)) * t71, -g(1) * t69 - g(2) * t70 - t43 * t68];
taug = t1(:);
