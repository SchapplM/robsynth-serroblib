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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:26
% EndTime: 2020-01-03 11:41:28
% DurationCPUTime: 0.43s
% Computational Cost: add. (233->86), mult. (307->108), div. (0->0), fcn. (285->10), ass. (0->54)
t67 = m(5) + m(6);
t74 = m(3) + m(4) + t67;
t68 = m(5) * pkin(3);
t72 = -mrSges(4,1) - t68;
t71 = mrSges(2,2) - mrSges(3,3);
t32 = qJ(3) + pkin(9);
t25 = cos(t32);
t38 = cos(qJ(3));
t29 = t38 * pkin(3);
t19 = pkin(4) * t25 + t29;
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t35 = -qJ(4) - pkin(6);
t69 = -mrSges(2,1) + (-mrSges(3,1) - m(4) * pkin(2) - m(5) * (t29 + pkin(2)) - m(6) * (pkin(2) + t19)) * t34 + (mrSges(3,2) - m(4) * pkin(6) - mrSges(4,3) + m(5) * t35 - mrSges(5,3) - m(6) * (pkin(7) - t35) - mrSges(6,3)) * t33;
t26 = qJ(5) + t32;
t22 = cos(t26);
t39 = cos(qJ(1));
t53 = t39 * t22;
t21 = sin(t26);
t37 = sin(qJ(1));
t61 = t37 * t21;
t5 = -t34 * t61 - t53;
t54 = t39 * t21;
t60 = t37 * t22;
t6 = t34 * t60 - t54;
t66 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t7 = t34 * t54 - t60;
t8 = t34 * t53 + t61;
t65 = t7 * mrSges(6,1) + t8 * mrSges(6,2);
t64 = g(1) * t33;
t36 = sin(qJ(3));
t63 = t36 * pkin(3);
t24 = sin(t32);
t18 = pkin(4) * t24 + t63;
t62 = t37 * t18;
t59 = t37 * t24;
t58 = t37 * t25;
t57 = t37 * t36;
t56 = t37 * t38;
t55 = t39 * t18;
t52 = t39 * t24;
t51 = t39 * t25;
t50 = t39 * t36;
t49 = t39 * t38;
t43 = -mrSges(6,1) * t21 - mrSges(6,2) * t22;
t15 = t34 * t50 - t56;
t13 = -t34 * t57 - t49;
t16 = t34 * t49 + t57;
t14 = t34 * t56 - t50;
t12 = t34 * t51 + t59;
t11 = t34 * t52 - t58;
t10 = t34 * t58 - t52;
t9 = -t34 * t59 - t51;
t1 = [(t50 * t68 + m(6) * t55 - t14 * mrSges(4,1) - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t13 * mrSges(4,2) - t9 * mrSges(5,2) - t5 * mrSges(6,2) - t74 * (t37 * pkin(1) - t39 * qJ(2)) - t71 * t39 + t69 * t37) * g(3) + (-t57 * t68 - m(6) * t62 - t16 * mrSges(4,1) - t12 * mrSges(5,1) - t8 * mrSges(6,1) + t15 * mrSges(4,2) + t11 * mrSges(5,2) + t7 * mrSges(6,2) - t74 * (t39 * pkin(1) + t37 * qJ(2)) + t71 * t37 + t69 * t39) * g(2), (g(2) * t39 + g(3) * t37) * t74, (m(5) * t63 + m(6) * t18 + mrSges(4,1) * t36 + mrSges(5,1) * t24 + mrSges(4,2) * t38 + mrSges(5,2) * t25 - t43) * t64 + (-t16 * mrSges(4,2) - t11 * mrSges(5,1) - t12 * mrSges(5,2) - m(6) * (-t37 * t19 + t34 * t55) - t65 + t72 * t15) * g(3) + (t14 * mrSges(4,2) - t9 * mrSges(5,1) + t10 * mrSges(5,2) - m(6) * (-t39 * t19 - t34 * t62) - t66 + t72 * t13) * g(2), (t34 * g(1) + t33 * (-g(2) * t37 + g(3) * t39)) * t67, -g(2) * t66 - g(3) * t65 - t43 * t64];
taug = t1(:);
