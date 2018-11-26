% Calculate Gravitation load on the joints for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2018-11-23 16:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:51:47
% EndTime: 2018-11-23 16:51:48
% DurationCPUTime: 1.06s
% Computational Cost: add. (1034->127), mult. (1220->155), div. (0->0), fcn. (1149->14), ass. (0->71)
t112 = m(6) + m(7);
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t106 = m(7) * pkin(5) + t49 * mrSges(7,1) - t45 * mrSges(7,2) + mrSges(6,1);
t102 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t98 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3) + t112 * (qJ(3) - pkin(9));
t111 = -t45 * mrSges(7,1) - t49 * mrSges(7,2) + t98;
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t86 = pkin(6) + qJ(2);
t67 = cos(t86) / 0.2e1;
t87 = pkin(6) - qJ(2);
t73 = cos(t87);
t56 = t73 / 0.2e1 + t67;
t97 = cos(qJ(1));
t23 = t47 * t48 - t97 * t56;
t71 = sin(t86);
t65 = t71 / 0.2e1;
t72 = sin(t87);
t62 = t65 - t72 / 0.2e1;
t51 = cos(qJ(2));
t93 = t48 * t51;
t24 = t97 * t62 + t93;
t46 = sin(qJ(5));
t50 = cos(qJ(5));
t43 = sin(pkin(6));
t79 = t43 * t97;
t4 = t24 * t50 + t46 * t79;
t110 = t23 * t49 + t4 * t45;
t109 = t23 * t45 - t4 * t49;
t88 = m(5) + t112;
t105 = mrSges(3,1) + mrSges(4,1) + mrSges(5,1);
t103 = -mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t101 = -t102 * t46 + t106 * t50 + t105;
t94 = t43 * t48;
t91 = t97 * pkin(1) + pkin(8) * t94;
t26 = t97 * t47 + t48 * t56;
t90 = qJ(3) * t26;
t89 = t23 * qJ(3);
t78 = t97 * t51;
t27 = -t48 * t62 + t78;
t83 = t27 * pkin(2) + t91;
t77 = -t48 * pkin(1) + pkin(8) * t79;
t11 = t23 * pkin(2);
t66 = t72 / 0.2e1;
t55 = t66 - t71 / 0.2e1;
t25 = -t97 * t55 + t93;
t76 = qJ(3) * t25 - t11;
t17 = t26 * pkin(2);
t28 = t48 * t55 + t78;
t75 = qJ(3) * t28 - t17;
t34 = t65 + t66;
t33 = t34 * pkin(2);
t35 = t67 - t73 / 0.2e1;
t74 = -qJ(3) * t35 + t33;
t70 = -t24 * pkin(2) + t77;
t60 = -t24 * t46 + t50 * t79;
t59 = t27 * pkin(3) - qJ(4) * t94 + t83;
t58 = t27 * pkin(4) + t59;
t54 = -t24 * pkin(3) - qJ(4) * t79 + t70;
t53 = -t24 * pkin(4) + t54;
t44 = cos(pkin(6));
t32 = t34 * pkin(3);
t22 = -t35 * t50 - t44 * t46;
t16 = t26 * pkin(3);
t10 = t23 * pkin(3);
t8 = t27 * t50 - t46 * t94;
t7 = t27 * t46 + t50 * t94;
t2 = -t26 * t45 + t49 * t8;
t1 = -t26 * t49 - t45 * t8;
t3 = [(-t97 * mrSges(2,1) + t48 * mrSges(2,2) - m(3) * t91 - m(4) * (t83 + t90) - m(5) * (t59 + t90) - m(6) * t58 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t58) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t102 * t7 + t103 * t94 - t105 * t27 - t98 * t26) * g(2) + (t48 * mrSges(2,1) + t97 * mrSges(2,2) - m(3) * t77 - m(4) * (t70 - t89) - m(5) * (t54 - t89) - m(6) * t53 + t4 * mrSges(6,1) - m(7) * (-pkin(5) * t4 + t53) - t109 * mrSges(7,1) - t110 * mrSges(7,2) + t102 * t60 + t103 * t79 + t105 * t24 + t98 * t23) * g(1) (-m(4) * t74 - m(5) * (t32 + t74) - t112 * (t34 * pkin(4) + t32 + t33) + t111 * t35 - t101 * t34) * g(3) + (-m(4) * t76 - m(5) * (-t10 + t76) - t112 * (-t23 * pkin(4) - t10 - t11) - t111 * t25 + t101 * t23) * g(2) + (-m(4) * t75 - m(5) * (-t16 + t75) - t112 * (-t26 * pkin(4) - t16 - t17) - t111 * t28 + t101 * t26) * g(1) (-g(1) * t26 - g(2) * t23 + g(3) * t34) * (m(4) + t88) (t44 * g(3) + (t48 * g(1) - t97 * g(2)) * t43) * t88 (t102 * t22 - t106 * (t35 * t46 - t44 * t50)) * g(3) + (t102 * t4 - t106 * t60) * g(2) + (t102 * t8 + t106 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t110 * mrSges(7,1) + t109 * mrSges(7,2)) - g(3) * ((-t22 * t45 + t34 * t49) * mrSges(7,1) + (-t22 * t49 - t34 * t45) * mrSges(7,2))];
taug  = t3(:);
