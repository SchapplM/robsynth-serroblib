% Calculate Gravitation load on the joints for
% S6RPRPRR3
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
% Datum: 2018-11-23 16:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:04:03
% EndTime: 2018-11-23 16:04:04
% DurationCPUTime: 0.56s
% Computational Cost: add. (456->101), mult. (412->116), div. (0->0), fcn. (367->12), ass. (0->54)
t79 = mrSges(6,3) + mrSges(7,3);
t32 = cos(pkin(11));
t20 = t32 * pkin(4) + pkin(3);
t29 = pkin(11) + qJ(5);
t23 = cos(t29);
t13 = pkin(5) * t23 + t20;
t25 = qJ(6) + t29;
t18 = sin(t25);
t19 = cos(t25);
t21 = sin(t29);
t78 = -m(6) * t20 - m(7) * t13 - t23 * mrSges(6,1) - t19 * mrSges(7,1) + t21 * mrSges(6,2) + t18 * mrSges(7,2);
t33 = -pkin(8) - qJ(4);
t28 = -pkin(9) + t33;
t77 = m(6) * t33 + m(7) * t28 - t79;
t56 = m(5) + m(6) + m(7);
t73 = -m(4) - t56;
t30 = qJ(1) + pkin(10);
t22 = sin(t30);
t24 = cos(t30);
t76 = g(1) * t24 + g(2) * t22;
t36 = cos(qJ(3));
t31 = sin(pkin(11));
t44 = m(5) * pkin(3) + t32 * mrSges(5,1) - t31 * mrSges(5,2);
t75 = t44 * t36;
t74 = m(7) * pkin(5) + mrSges(6,1);
t34 = sin(qJ(3));
t50 = t36 * mrSges(4,1) - t34 * mrSges(4,2);
t72 = t79 * t34 + mrSges(3,1) + t50;
t62 = t31 * pkin(4);
t66 = pkin(5) * t21;
t70 = -m(6) * t62 - m(7) * (t62 + t66) + mrSges(3,2) - mrSges(4,3) - t31 * mrSges(5,1) - t32 * mrSges(5,2);
t60 = t22 * t36;
t5 = t18 * t60 + t24 * t19;
t6 = t24 * t18 - t19 * t60;
t68 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t59 = t24 * t36;
t7 = -t18 * t59 + t22 * t19;
t8 = t22 * t18 + t19 * t59;
t67 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t63 = g(3) * t34;
t35 = sin(qJ(1));
t61 = t35 * pkin(1);
t37 = cos(qJ(1));
t27 = t37 * pkin(1);
t51 = m(5) * qJ(4) + mrSges(5,3);
t47 = -mrSges(7,1) * t18 - mrSges(7,2) * t19;
t46 = t36 * t13 - t34 * t28;
t45 = t36 * t20 - t34 * t33;
t11 = -t21 * t59 + t22 * t23;
t9 = t21 * t60 + t24 * t23;
t41 = t51 * t34 + t75;
t12 = t22 * t21 + t23 * t59;
t10 = t24 * t21 - t23 * t60;
t1 = [(-m(3) * t27 - t37 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t35 * mrSges(2,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + t73 * (t24 * pkin(2) + t22 * pkin(7) + t27) + t70 * t22 + (-m(6) * t45 - m(7) * t46 - t41 - t72) * t24) * g(2) + (m(3) * t61 + t35 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t37 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t73 * (t24 * pkin(7) - t61) + t70 * t24 + (m(4) * pkin(2) - m(5) * (-t34 * qJ(4) - pkin(2)) + t34 * mrSges(5,3) + t75 - m(6) * (-pkin(2) - t45) - m(7) * (-pkin(2) - t46) + t72) * t22) * g(1) (-m(3) + t73) * g(3) (-t41 - t50) * g(3) + (t78 * g(3) + t76 * (mrSges(4,2) - t51 + t77)) * t36 + (t77 * g(3) + t76 * (mrSges(4,1) + t44 - t78)) * t34 (t36 * g(3) - t76 * t34) * t56 (m(7) * t66 + mrSges(6,1) * t21 + mrSges(6,2) * t23 - t47) * t63 + (-t10 * mrSges(6,2) + t74 * t9 - t68) * g(2) + (t12 * mrSges(6,2) - t74 * t11 - t67) * g(1), -g(1) * t67 - g(2) * t68 - t47 * t63];
taug  = t1(:);
