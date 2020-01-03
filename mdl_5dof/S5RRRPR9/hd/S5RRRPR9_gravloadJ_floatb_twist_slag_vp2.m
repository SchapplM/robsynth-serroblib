% Calculate Gravitation load on the joints for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:15
% EndTime: 2019-12-31 21:22:17
% DurationCPUTime: 0.63s
% Computational Cost: add. (315->107), mult. (424->129), div. (0->0), fcn. (389->10), ass. (0->61)
t86 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t32 = qJ(3) + pkin(9);
t25 = cos(t32);
t37 = cos(qJ(3));
t28 = t37 * pkin(3);
t19 = pkin(4) * t25 + t28;
t17 = pkin(2) + t19;
t26 = qJ(5) + t32;
t21 = sin(t26);
t22 = cos(t26);
t23 = t28 + pkin(2);
t24 = sin(t32);
t34 = sin(qJ(3));
t85 = -m(4) * pkin(2) - m(5) * t23 - m(6) * t17 - t37 * mrSges(4,1) - t25 * mrSges(5,1) - t22 * mrSges(6,1) + t34 * mrSges(4,2) + t24 * mrSges(5,2) + t21 * mrSges(6,2);
t33 = -qJ(4) - pkin(7);
t31 = -pkin(8) + t33;
t84 = -m(4) * pkin(7) + m(5) * t33 + m(6) * t31 - t86;
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t83 = g(1) * t39 + g(2) * t36;
t75 = m(5) * pkin(3);
t82 = m(5) + m(6);
t69 = t34 * pkin(3);
t18 = pkin(4) * t24 + t69;
t81 = m(6) * t18;
t80 = mrSges(4,1) + t75;
t79 = mrSges(2,2) - mrSges(3,3);
t78 = -m(3) - m(4) - t82;
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t50 = t38 * mrSges(3,1) - t35 * mrSges(3,2);
t76 = t86 * t35 + mrSges(2,1) + t50;
t61 = t39 * t22;
t64 = t36 * t38;
t5 = t21 * t64 + t61;
t62 = t39 * t21;
t6 = -t22 * t64 + t62;
t74 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = t36 * t22 - t38 * t62;
t8 = t36 * t21 + t38 * t61;
t73 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t70 = g(3) * t35;
t65 = t36 * t34;
t63 = t39 * t18;
t60 = t39 * t24;
t59 = t39 * t25;
t58 = t39 * t34;
t57 = t39 * t37;
t52 = t38 * pkin(2) + t35 * pkin(7);
t48 = -mrSges(6,1) * t21 - mrSges(6,2) * t22;
t47 = t38 * t17 - t35 * t31;
t46 = t38 * t23 - t35 * t33;
t15 = t36 * t37 - t38 * t58;
t13 = t34 * t64 + t57;
t16 = t38 * t57 + t65;
t14 = -t37 * t64 + t58;
t12 = t36 * t24 + t38 * t59;
t11 = t36 * t25 - t38 * t60;
t10 = -t25 * t64 + t60;
t9 = t24 * t64 + t59;
t1 = [(-t65 * t75 - t16 * mrSges(4,1) - t12 * mrSges(5,1) - t8 * mrSges(6,1) - t15 * mrSges(4,2) - t11 * mrSges(5,2) - t7 * mrSges(6,2) + t78 * (t39 * pkin(1) + t36 * pkin(6)) + (t79 - t81) * t36 + (-m(4) * t52 - m(5) * t46 - m(6) * t47 - t76) * t39) * g(2) + (-t58 * t75 - m(6) * t63 - t14 * mrSges(4,1) - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t13 * mrSges(4,2) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t52) - m(5) * (-pkin(1) - t46) - m(6) * (-pkin(1) - t47) + t76) * t36 + (t78 * pkin(6) + t79) * t39) * g(1), -g(3) * t50 + (t85 * g(3) + t83 * (mrSges(3,2) + t84)) * t38 + (t84 * g(3) + t83 * (mrSges(3,1) - t85)) * t35, (m(5) * t69 + mrSges(4,1) * t34 + mrSges(5,1) * t24 + mrSges(4,2) * t37 + mrSges(5,2) * t25 - t48 + t81) * t70 + (-t14 * mrSges(4,2) + t9 * mrSges(5,1) - t10 * mrSges(5,2) - m(6) * (-t18 * t64 - t39 * t19) - t74 + t80 * t13) * g(2) + (t16 * mrSges(4,2) - t11 * mrSges(5,1) + t12 * mrSges(5,2) - m(6) * (t36 * t19 - t38 * t63) - t73 - t80 * t15) * g(1), (t38 * g(3) - t35 * t83) * t82, -g(1) * t73 - g(2) * t74 - t48 * t70];
taug = t1(:);
