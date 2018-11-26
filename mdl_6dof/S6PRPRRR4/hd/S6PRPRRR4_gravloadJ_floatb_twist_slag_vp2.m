% Calculate Gravitation load on the joints for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:05:22
% EndTime: 2018-11-23 15:05:23
% DurationCPUTime: 0.82s
% Computational Cost: add. (1133->92), mult. (1125->123), div. (0->0), fcn. (1079->18), ass. (0->47)
t41 = qJ(5) + qJ(6);
t38 = sin(t41);
t39 = cos(t41);
t48 = sin(qJ(5));
t50 = cos(qJ(5));
t99 = mrSges(5,1) + m(7) * (pkin(5) * t50 + pkin(4)) + t39 * mrSges(7,1) - t38 * mrSges(7,2) + m(6) * pkin(4) + t50 * mrSges(6,1) - t48 * mrSges(6,2);
t89 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t40 = pkin(12) + qJ(4);
t36 = sin(t40);
t37 = cos(t40);
t45 = cos(pkin(12));
t88 = m(4) * pkin(2) + t45 * mrSges(4,1) - sin(pkin(12)) * mrSges(4,2) - t89 * t36 + mrSges(3,1) + t99 * t37;
t95 = -m(7) * pkin(5) - mrSges(6,1);
t92 = -m(5) - m(6) - m(7);
t87 = -m(4) * qJ(3) - mrSges(7,1) * t38 - mrSges(6,2) * t50 - mrSges(7,2) * t39 + t95 * t48 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t43 = sin(pkin(11));
t49 = sin(qJ(2));
t76 = pkin(6) + qJ(2);
t67 = cos(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t70 = cos(t77);
t55 = t70 / 0.2e1 + t67;
t78 = cos(pkin(11));
t17 = t43 * t49 - t55 * t78;
t66 = sin(t76) / 0.2e1;
t69 = sin(t77);
t28 = t66 - t69 / 0.2e1;
t51 = cos(qJ(2));
t61 = t28 * t78 + t43 * t51;
t44 = sin(pkin(6));
t71 = t44 * t78;
t8 = -t36 * t71 + t37 * t61;
t85 = (t17 * t39 - t38 * t8) * mrSges(7,1) + (-t17 * t38 - t39 * t8) * mrSges(7,2);
t60 = -t43 * t28 + t51 * t78;
t82 = t43 * t44;
t10 = t36 * t82 + t37 * t60;
t20 = t43 * t55 + t49 * t78;
t84 = (-t10 * t38 + t20 * t39) * mrSges(7,1) + (-t10 * t39 - t20 * t38) * mrSges(7,2);
t29 = t67 - t70 / 0.2e1;
t46 = cos(pkin(6));
t14 = -t29 * t37 + t36 * t46;
t27 = t66 + t69 / 0.2e1;
t83 = (-t14 * t38 - t27 * t39) * mrSges(7,1) + (-t14 * t39 + t27 * t38) * mrSges(7,2);
t75 = m(4) - t92;
t47 = -pkin(8) - qJ(3);
t34 = pkin(3) * t45 + pkin(2);
t1 = [(-m(2) - m(3) - t75) * g(3) (t92 * (t27 * t34 + t29 * t47) - t87 * t29 - t88 * t27) * g(3) + (t92 * (-t17 * t34 - t47 * t61) + t87 * t61 + t88 * t17) * g(2) + (t92 * (-t20 * t34 - t47 * t60) + t87 * t60 + t88 * t20) * g(1) (-g(1) * t20 - g(2) * t17 + g(3) * t27) * t75 (t89 * t14 - t99 * (t29 * t36 + t37 * t46)) * g(3) + (t89 * t8 - t99 * (-t36 * t61 - t37 * t71)) * g(2) + (-t99 * (-t36 * t60 + t37 * t82) + t89 * t10) * g(1) (-(-t14 * t50 + t27 * t48) * mrSges(6,2) - t83 + t95 * (-t14 * t48 - t27 * t50)) * g(3) + (-(-t17 * t48 - t50 * t8) * mrSges(6,2) - t85 + t95 * (t17 * t50 - t48 * t8)) * g(2) + (-(-t10 * t50 - t20 * t48) * mrSges(6,2) - t84 + t95 * (-t10 * t48 + t20 * t50)) * g(1), -g(1) * t84 - g(2) * t85 - g(3) * t83];
taug  = t1(:);
