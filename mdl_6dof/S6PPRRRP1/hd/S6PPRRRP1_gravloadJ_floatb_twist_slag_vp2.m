% Calculate Gravitation load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:21
% EndTime: 2019-03-08 18:52:23
% DurationCPUTime: 0.76s
% Computational Cost: add. (808->87), mult. (2225->139), div. (0->0), fcn. (2814->14), ass. (0->61)
t104 = mrSges(6,1) + mrSges(7,1);
t99 = -mrSges(6,2) - mrSges(7,2);
t44 = cos(qJ(5));
t103 = m(6) * pkin(4) + m(7) * (pkin(5) * t44 + pkin(4)) + mrSges(5,1);
t102 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t89 = m(7) * pkin(5);
t100 = mrSges(4,2) - mrSges(5,3);
t75 = sin(pkin(12));
t76 = sin(pkin(11));
t60 = t76 * t75;
t79 = cos(pkin(12));
t80 = cos(pkin(11));
t67 = t80 * t79;
t82 = cos(pkin(6));
t51 = -t82 * t67 + t60;
t77 = sin(pkin(7));
t78 = sin(pkin(6));
t64 = t78 * t77;
t81 = cos(pkin(7));
t98 = t51 * t81 + t80 * t64;
t61 = t76 * t79;
t65 = t80 * t75;
t52 = t82 * t61 + t65;
t63 = t78 * t76;
t97 = t52 * t81 - t77 * t63;
t96 = t79 * t81 * t78 + t82 * t77;
t95 = -m(5) - m(6) - m(7);
t41 = sin(qJ(5));
t94 = t104 * t44 + t99 * t41 + t103;
t93 = -t89 - t104;
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t91 = t102 * t42 + t103 * t45 + mrSges(4,1);
t90 = m(3) + m(4) - t95;
t88 = cos(qJ(3));
t33 = t82 * t65 + t61;
t43 = sin(qJ(3));
t20 = t33 * t88 - t98 * t43;
t87 = t20 * t41;
t34 = -t82 * t60 + t67;
t22 = t34 * t88 - t97 * t43;
t86 = t22 * t41;
t62 = t78 * t75;
t29 = t96 * t43 + t88 * t62;
t85 = t29 * t41;
t84 = t41 * t45;
t83 = t44 * t45;
t66 = t80 * t78;
t50 = -t79 * t64 + t82 * t81;
t47 = t52 * t77 + t81 * t63;
t46 = t51 * t77 - t81 * t66;
t28 = t43 * t62 - t96 * t88;
t24 = t29 * t45 + t50 * t42;
t23 = t29 * t42 - t50 * t45;
t21 = t34 * t43 + t97 * t88;
t19 = t33 * t43 + t98 * t88;
t12 = t22 * t45 + t47 * t42;
t11 = t22 * t42 - t47 * t45;
t10 = t20 * t45 + t46 * t42;
t9 = t20 * t42 - t46 * t45;
t1 = [(-m(2) - t90) * g(3) (-t63 * g(1) + t66 * g(2) - t82 * g(3)) * t90 (-t85 * t89 + t95 * (-t28 * pkin(3) + t29 * pkin(9)) + t100 * t29 - t104 * (-t28 * t83 + t85) + t99 * (t28 * t84 + t29 * t44) + t91 * t28) * g(3) + (-t87 * t89 + t95 * (-t19 * pkin(3) + t20 * pkin(9)) - t104 * (-t19 * t83 + t87) + t99 * (t19 * t84 + t20 * t44) + t100 * t20 + t91 * t19) * g(2) + (-t86 * t89 - t104 * (-t21 * t83 + t86) + t95 * (-t21 * pkin(3) + t22 * pkin(9)) + t99 * (t21 * t84 + t22 * t44) + t100 * t22 + t91 * t21) * g(1) (-t102 * t24 + t94 * t23) * g(3) + (-t10 * t102 + t94 * t9) * g(2) + (-t102 * t12 + t94 * t11) * g(1) (t99 * (-t24 * t44 - t28 * t41) + t93 * (-t24 * t41 + t28 * t44)) * g(3) + (t99 * (-t10 * t44 - t19 * t41) + t93 * (-t10 * t41 + t19 * t44)) * g(2) + (t99 * (-t12 * t44 - t21 * t41) + t93 * (-t12 * t41 + t21 * t44)) * g(1) (-g(1) * t11 - g(2) * t9 - g(3) * t23) * m(7)];
taug  = t1(:);
