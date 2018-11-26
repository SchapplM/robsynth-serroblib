% Calculate Gravitation load on the joints for
% S6PRRPRR5
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
% Datum: 2018-11-23 15:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:17:24
% EndTime: 2018-11-23 15:17:25
% DurationCPUTime: 0.92s
% Computational Cost: add. (1229->98), mult. (1339->125), div. (0->0), fcn. (1311->18), ass. (0->57)
t42 = pkin(12) + qJ(5);
t39 = qJ(6) + t42;
t34 = sin(t39);
t35 = cos(t39);
t44 = cos(pkin(12));
t36 = t44 * pkin(4) + pkin(3);
t37 = sin(t42);
t38 = cos(t42);
t43 = sin(pkin(12));
t100 = mrSges(4,1) + m(7) * (pkin(5) * t38 + t36) + t35 * mrSges(7,1) - t34 * mrSges(7,2) + m(6) * t36 + t38 * mrSges(6,1) - t37 * mrSges(6,2) + m(5) * pkin(3) + t44 * mrSges(5,1) - t43 * mrSges(5,2);
t45 = -pkin(9) - qJ(4);
t99 = mrSges(4,2) - m(5) * qJ(4) - mrSges(5,3) + m(7) * (-pkin(10) + t45) - mrSges(7,3) + m(6) * t45 - mrSges(6,3);
t104 = m(5) + m(6) + m(7);
t101 = m(4) + t104;
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t108 = pkin(2) * t101 + t100 * t48 - t99 * t46 + mrSges(3,1);
t82 = pkin(4) * t43 + pkin(8);
t107 = m(6) * t82 + t37 * mrSges(6,1) + t38 * mrSges(6,2) + m(7) * (pkin(5) * t37 + t82) + t34 * mrSges(7,1) + t35 * mrSges(7,2) + t43 * mrSges(5,1) + t44 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3) - (-m(4) - m(5)) * pkin(8);
t105 = -m(7) * pkin(5) - mrSges(6,1);
t47 = sin(qJ(2));
t83 = pkin(6) + qJ(2);
t72 = cos(t83) / 0.2e1;
t84 = pkin(6) - qJ(2);
t76 = cos(t84);
t54 = t76 / 0.2e1 + t72;
t85 = sin(pkin(11));
t87 = cos(pkin(11));
t13 = t85 * t47 - t87 * t54;
t74 = sin(t83);
t70 = t74 / 0.2e1;
t75 = sin(t84);
t61 = t70 - t75 / 0.2e1;
t49 = cos(qJ(2));
t77 = t85 * t49;
t14 = t87 * t61 + t77;
t86 = sin(pkin(6));
t64 = t87 * t86;
t8 = t14 * t48 - t46 * t64;
t93 = (t13 * t35 - t34 * t8) * mrSges(7,1) + (-t13 * t34 - t35 * t8) * mrSges(7,2);
t78 = t87 * t49;
t17 = -t85 * t61 + t78;
t63 = t86 * t85;
t10 = t17 * t48 + t46 * t63;
t16 = t87 * t47 + t85 * t54;
t92 = (-t10 * t34 + t16 * t35) * mrSges(7,1) + (-t10 * t35 - t16 * t34) * mrSges(7,2);
t26 = t72 - t76 / 0.2e1;
t88 = cos(pkin(6));
t20 = -t26 * t48 + t88 * t46;
t71 = t75 / 0.2e1;
t25 = t70 + t71;
t91 = (-t20 * t34 - t25 * t35) * mrSges(7,1) + (-t20 * t35 + t25 * t34) * mrSges(7,2);
t62 = t71 - t74 / 0.2e1;
t19 = -t26 * t46 - t88 * t48;
t9 = t17 * t46 - t48 * t63;
t7 = t14 * t46 + t48 * t64;
t1 = [(-m(2) - m(3) - t101) * g(3) (t107 * t26 - t108 * t25) * g(3) + (-t107 * (-t87 * t62 + t77) + t108 * t13) * g(2) + (-t107 * (t85 * t62 + t78) + t108 * t16) * g(1) (t100 * t19 + t99 * t20) * g(3) + (t100 * t7 + t99 * t8) * g(2) + (t99 * t10 + t100 * t9) * g(1), t104 * (-g(1) * t9 - g(2) * t7 - g(3) * t19) (-(-t20 * t38 + t25 * t37) * mrSges(6,2) - t91 + t105 * (-t20 * t37 - t25 * t38)) * g(3) + (-(-t13 * t37 - t38 * t8) * mrSges(6,2) - t93 + t105 * (t13 * t38 - t37 * t8)) * g(2) + (-(-t10 * t38 - t16 * t37) * mrSges(6,2) - t92 + t105 * (-t10 * t37 + t16 * t38)) * g(1), -g(1) * t92 - g(2) * t93 - g(3) * t91];
taug  = t1(:);
