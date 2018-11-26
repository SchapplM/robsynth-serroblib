% Calculate Gravitation load on the joints for
% S6RPRPRR1
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

function taug = S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:00
% EndTime: 2018-11-23 16:03:00
% DurationCPUTime: 0.48s
% Computational Cost: add. (422->88), mult. (317->94), div. (0->0), fcn. (259->12), ass. (0->50)
t82 = -m(6) - m(7);
t31 = qJ(3) + pkin(11);
t27 = qJ(5) + t31;
t20 = sin(t27);
t21 = cos(t27);
t37 = cos(qJ(6));
t61 = t37 * mrSges(7,1);
t40 = mrSges(6,2) * t21 + (m(7) * pkin(5) + mrSges(6,1) + t61) * t20;
t34 = sin(qJ(6));
t62 = t34 * mrSges(7,2);
t90 = -t20 * t62 + t21 * (-m(7) * pkin(9) - mrSges(7,3));
t80 = t21 * pkin(5) + t20 * pkin(9);
t89 = m(7) * t80;
t88 = -t21 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t20;
t23 = sin(t31);
t25 = cos(t31);
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t72 = t35 * pkin(3);
t87 = m(5) * t72 + mrSges(4,1) * t35 + mrSges(5,1) * t23 + mrSges(4,2) * t38 + mrSges(5,2) * t25 + t82 * (-pkin(4) * t23 - t72) + t40;
t86 = -t38 * mrSges(4,1) - t25 * mrSges(5,1) + t35 * mrSges(4,2) + t23 * mrSges(5,2);
t83 = -m(3) - m(4);
t79 = -(t61 - t62) * t21 + t88;
t78 = m(5) - t83;
t28 = t38 * pkin(3);
t76 = mrSges(3,1) + m(5) * (t28 + pkin(2)) + m(4) * pkin(2) - t86 - t88;
t33 = -qJ(4) - pkin(7);
t75 = -m(4) * pkin(7) + m(5) * t33 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t36 = sin(qJ(1));
t71 = t36 * pkin(1);
t39 = cos(qJ(1));
t29 = t39 * pkin(1);
t32 = qJ(1) + pkin(10);
t24 = sin(t32);
t66 = t24 * t34;
t65 = t24 * t37;
t26 = cos(t32);
t64 = t26 * t34;
t63 = t26 * t37;
t59 = pkin(4) * t25 + t28;
t58 = m(5) - t82;
t55 = t90 * t24;
t54 = t90 * t26;
t30 = -pkin(8) + t33;
t12 = pkin(2) + t59;
t4 = t21 * t63 + t66;
t3 = -t21 * t64 + t65;
t2 = -t21 * t65 + t64;
t1 = t21 * t66 + t63;
t5 = [(-t39 * mrSges(2,1) - t4 * mrSges(7,1) + t36 * mrSges(2,2) - t3 * mrSges(7,2) + t82 * (t26 * t12 - t24 * t30 + t29) - t78 * t29 + t75 * t24 + (-t76 - t89) * t26) * g(2) + (t36 * mrSges(2,1) - t2 * mrSges(7,1) + t39 * mrSges(2,2) - t1 * mrSges(7,2) + t82 * (-t26 * t30 - t71) + t78 * t71 + t75 * t26 + (m(6) * t12 - m(7) * (-t12 - t80) + t76) * t24) * g(1) (-t58 + t83) * g(3) (-m(5) * t28 - m(6) * t59 - m(7) * (t59 + t80) + t79 + t86) * g(3) + (t87 * t24 + t55) * g(2) + (t87 * t26 + t54) * g(1) (-g(1) * t24 + g(2) * t26) * t58 (t79 - t89) * g(3) + (t40 * t24 + t55) * g(2) + (t40 * t26 + t54) * g(1), -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t34 - mrSges(7,2) * t37) * t20];
taug  = t5(:);
