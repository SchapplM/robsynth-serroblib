% Calculate Gravitation load on the joints for
% S6RPRPRR6
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
% Datum: 2018-11-23 16:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:05:42
% EndTime: 2018-11-23 16:05:43
% DurationCPUTime: 0.64s
% Computational Cost: add. (448->109), mult. (435->117), div. (0->0), fcn. (385->12), ass. (0->58)
t82 = mrSges(6,3) + mrSges(7,3);
t31 = cos(pkin(11));
t18 = t31 * pkin(4) + pkin(3);
t27 = pkin(11) + qJ(5);
t22 = cos(t27);
t13 = pkin(5) * t22 + t18;
t24 = qJ(6) + t27;
t16 = sin(t24);
t17 = cos(t24);
t20 = sin(t27);
t81 = -m(6) * t18 - m(7) * t13 - t22 * mrSges(6,1) - t17 * mrSges(7,1) + t20 * mrSges(6,2) + t16 * mrSges(7,2);
t33 = -pkin(8) - qJ(4);
t26 = -pkin(9) + t33;
t80 = m(6) * t33 + m(7) * t26 - t82;
t35 = sin(qJ(1));
t36 = cos(qJ(1));
t79 = g(1) * t36 + g(2) * t35;
t78 = -m(6) - m(7);
t28 = pkin(10) + qJ(3);
t23 = cos(t28);
t29 = sin(pkin(11));
t44 = m(5) * pkin(3) + t31 * mrSges(5,1) - t29 * mrSges(5,2);
t77 = t44 * t23;
t76 = m(7) * pkin(5) + mrSges(6,1);
t34 = -pkin(7) - qJ(2);
t75 = -m(3) * qJ(2) + m(5) * t34 - t29 * mrSges(5,1) - t31 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t21 = sin(t28);
t32 = cos(pkin(10));
t50 = t23 * mrSges(4,1) - t21 * mrSges(4,2);
t73 = m(3) * pkin(1) + t32 * mrSges(3,1) - sin(pkin(10)) * mrSges(3,2) + mrSges(2,1) + t50 + t82 * t21;
t61 = t35 * t16;
t5 = t17 * t36 + t23 * t61;
t60 = t35 * t17;
t6 = t16 * t36 - t23 * t60;
t71 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t62 = t23 * t36;
t7 = -t16 * t62 + t60;
t8 = t17 * t62 + t61;
t70 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t69 = pkin(4) * t29;
t68 = pkin(5) * t20;
t65 = g(3) * t21;
t59 = t35 * t20;
t58 = t35 * t22;
t57 = m(5) - t78;
t53 = m(5) * qJ(4) + mrSges(5,3);
t48 = -mrSges(7,1) * t16 - mrSges(7,2) * t17;
t47 = t13 * t23 - t21 * t26;
t46 = t18 * t23 - t21 * t33;
t11 = -t20 * t62 + t58;
t9 = t22 * t36 + t23 * t59;
t40 = t53 * t21 + t77;
t19 = pkin(2) * t32 + pkin(1);
t15 = t36 * t19;
t14 = t68 + t69;
t12 = t22 * t62 + t59;
t10 = t20 * t36 - t23 * t58;
t1 = [(-m(5) * t15 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + (-m(4) + t78) * (-t35 * t34 + t15) + (-m(6) * t69 - m(7) * t14 + t75) * t35 + (-m(6) * t46 - m(7) * t47 - t40 - t73) * t36) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (m(4) * t34 - m(6) * (-t34 + t69) - m(7) * (t14 - t34) + t75) * t36 + (m(4) * t19 - m(5) * (-qJ(4) * t21 - t19) + t21 * mrSges(5,3) + t77 - m(6) * (-t19 - t46) - m(7) * (-t19 - t47) + t73) * t35) * g(1) (-g(1) * t35 + g(2) * t36) * (m(3) + m(4) + t57) (-t40 - t50) * g(3) + (t81 * g(3) + t79 * (mrSges(4,2) - t53 + t80)) * t23 + (t80 * g(3) + t79 * (mrSges(4,1) + t44 - t81)) * t21 (t23 * g(3) - t79 * t21) * t57 (m(7) * t68 + mrSges(6,1) * t20 + mrSges(6,2) * t22 - t48) * t65 + (-mrSges(6,2) * t10 + t76 * t9 - t71) * g(2) + (mrSges(6,2) * t12 - t76 * t11 - t70) * g(1), -g(1) * t70 - g(2) * t71 - t48 * t65];
taug  = t1(:);
