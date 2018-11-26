% Calculate Gravitation load on the joints for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:18
% EndTime: 2018-11-23 16:03:19
% DurationCPUTime: 0.59s
% Computational Cost: add. (434->95), mult. (366->107), div. (0->0), fcn. (320->12), ass. (0->56)
t31 = cos(qJ(5));
t14 = pkin(5) * t31 + pkin(4);
t26 = qJ(5) + qJ(6);
t20 = sin(t26);
t21 = cos(t26);
t28 = sin(qJ(5));
t78 = -m(6) * pkin(4) - m(7) * t14 - mrSges(6,1) * t31 - mrSges(7,1) * t21 + mrSges(6,2) * t28 + mrSges(7,2) * t20;
t49 = m(5) + m(6) + m(7);
t77 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t24 = qJ(3) + pkin(11);
t16 = sin(t24);
t18 = cos(t24);
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t76 = -t32 * mrSges(4,1) - t18 * mrSges(5,1) + t29 * mrSges(4,2) - t16 * t77;
t69 = m(7) * pkin(5);
t75 = m(3) + m(4);
t74 = mrSges(6,1) + t69;
t71 = m(4) * pkin(2) + mrSges(3,1) - t76;
t70 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t25 = qJ(1) + pkin(10);
t19 = cos(t25);
t52 = t19 * t21;
t17 = sin(t25);
t57 = t17 * t20;
t5 = t18 * t57 + t52;
t53 = t19 * t20;
t56 = t17 * t21;
t6 = -t18 * t56 + t53;
t68 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t7 = -t18 * t53 + t56;
t8 = t18 * t52 + t57;
t67 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t63 = g(3) * t16;
t62 = t16 * pkin(8);
t30 = sin(qJ(1));
t61 = t30 * pkin(1);
t22 = t32 * pkin(3);
t33 = cos(qJ(1));
t23 = t33 * pkin(1);
t34 = -pkin(9) - pkin(8);
t58 = t16 * t34;
t55 = t17 * t28;
t54 = t17 * t31;
t51 = t19 * t28;
t50 = t19 * t31;
t47 = pkin(4) * t18 + t62;
t42 = -mrSges(7,1) * t20 - mrSges(7,2) * t21;
t41 = t14 * t18 - t58;
t11 = -t18 * t51 + t54;
t9 = t18 * t55 + t50;
t27 = -qJ(4) - pkin(7);
t15 = t22 + pkin(2);
t12 = t18 * t50 + t55;
t10 = -t18 * t54 + t51;
t1 = [(-t55 * t69 - t33 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t30 * mrSges(2,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t49 * (t15 * t19 - t17 * t27 + t23) - t75 * t23 + t70 * t17 + (-m(6) * t47 - m(7) * t41 - t71) * t19) * g(2) + (-t51 * t69 + t30 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t33 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t75 * t61 - t49 * (-t19 * t27 - t61) + t70 * t19 + (m(5) * t15 - m(6) * (-t15 - t47) - m(7) * (-t15 - t41) + t71) * t17) * g(1) (-t49 - t75) * g(3) (-m(5) * t22 - m(6) * (t22 + t62) - m(7) * (t22 - t58) + t78 * t18 + t76) * g(3) + (g(1) * t19 + g(2) * t17) * (mrSges(4,2) * t32 + (-m(6) * pkin(8) + m(7) * t34 - t77) * t18 + (mrSges(5,1) - t78) * t16 + (pkin(3) * t49 + mrSges(4,1)) * t29) (-g(1) * t17 + g(2) * t19) * t49 (mrSges(6,2) * t31 + t28 * t74 - t42) * t63 + (-t10 * mrSges(6,2) + t74 * t9 - t68) * g(2) + (t12 * mrSges(6,2) - t11 * t74 - t67) * g(1), -g(1) * t67 - g(2) * t68 - t42 * t63];
taug  = t1(:);
