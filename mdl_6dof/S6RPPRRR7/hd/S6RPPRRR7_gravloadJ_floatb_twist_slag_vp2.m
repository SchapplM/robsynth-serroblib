% Calculate Gravitation load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 15:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:51:10
% EndTime: 2018-11-23 15:51:10
% DurationCPUTime: 0.46s
% Computational Cost: add. (308->90), mult. (317->103), div. (0->0), fcn. (258->10), ass. (0->53)
t28 = pkin(10) + qJ(4);
t23 = qJ(5) + t28;
t19 = sin(t23);
t20 = cos(t23);
t83 = t19 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t20;
t82 = -m(4) - m(5);
t81 = -m(6) - m(7);
t49 = -t19 * pkin(5) + t20 * pkin(9);
t80 = m(7) * t49;
t35 = cos(qJ(1));
t32 = sin(qJ(6));
t64 = mrSges(7,2) * t32;
t53 = t20 * t64;
t65 = mrSges(6,2) * t19;
t79 = (-t53 - t65) * t35;
t33 = sin(qJ(1));
t34 = cos(qJ(6));
t59 = t34 * mrSges(7,1);
t52 = t20 * t59;
t62 = t33 * t20;
t63 = t19 * t33;
t78 = -mrSges(6,1) * t62 - mrSges(7,3) * t63 - (t52 - t53) * t33;
t77 = -(-t59 + t64) * t19 + t83;
t76 = mrSges(6,1) * t20 + t19 * mrSges(7,3) + t52;
t75 = -g(1) * t33 + g(2) * t35;
t74 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t29 = sin(pkin(10));
t21 = sin(t28);
t22 = cos(t28);
t44 = -t21 * mrSges(5,1) - t22 * mrSges(5,2);
t67 = t29 * pkin(3);
t73 = -m(5) * t67 + mrSges(2,2) - mrSges(3,3) + t44 - t29 * mrSges(4,1) - cos(pkin(10)) * mrSges(4,2) + t80 - t83;
t72 = pkin(4) * t21;
t71 = pkin(4) * t22;
t61 = t33 * t32;
t60 = t33 * t34;
t58 = t35 * t32;
t57 = t35 * t34;
t31 = -pkin(7) - qJ(3);
t56 = pkin(5) * t62 + pkin(9) * t63;
t55 = t35 * pkin(1) + t33 * qJ(2);
t54 = m(6) * t71;
t51 = t81 + t82;
t25 = t35 * qJ(2);
t50 = -t33 * pkin(1) + t25;
t47 = -pkin(5) * t20 - pkin(9) * t19;
t27 = -pkin(8) + t31;
t9 = t67 + t72;
t4 = t19 * t57 - t61;
t3 = t19 * t58 + t60;
t2 = t19 * t60 + t58;
t1 = -t19 * t61 + t57;
t5 = [(-t2 * mrSges(7,1) - t1 * mrSges(7,2) + t81 * (-t35 * t27 + t33 * t9 + t55) + (-m(3) + t82) * t55 + (-m(4) * qJ(3) + m(5) * t31 - t74) * t35 + t73 * t33) * g(2) + (-m(3) * t50 - t4 * mrSges(7,1) + t3 * mrSges(7,2) + t81 * (t33 * t27 + t35 * t9 + t50) + t82 * t25 + (-m(4) * (-pkin(1) - qJ(3)) - m(5) * (-pkin(1) + t31) + t74) * t33 + t73 * t35) * g(1), t75 * (m(3) - t51) (g(1) * t35 + g(2) * t33) * t51, t75 * (mrSges(5,1) * t22 - mrSges(5,2) * t21) + ((t54 - m(7) * (t47 - t71) + t76) * t35 + t79) * g(2) + (-(t54 - t65) * t33 - m(7) * (t33 * t71 + t56) + t78) * g(1) + (-t44 + m(6) * t72 - m(7) * (t49 - t72) + t77) * g(3) (t77 - t80) * g(3) + ((-m(7) * t47 + t76) * t35 + t79) * g(2) + (-m(7) * t56 + mrSges(6,2) * t63 + t78) * g(1), -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t32 - mrSges(7,2) * t34) * t20];
taug  = t5(:);
