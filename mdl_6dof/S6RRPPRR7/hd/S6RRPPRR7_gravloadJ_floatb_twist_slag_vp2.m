% Calculate Gravitation load on the joints for
% S6RRPPRR7
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
% Datum: 2018-11-23 16:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:53:16
% EndTime: 2018-11-23 16:53:17
% DurationCPUTime: 0.81s
% Computational Cost: add. (1034->119), mult. (1220->149), div. (0->0), fcn. (1149->14), ass. (0->68)
t96 = m(6) + m(7);
t100 = -t96 * pkin(9) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t47 = sin(qJ(6));
t51 = cos(qJ(6));
t97 = t47 * mrSges(7,1) + t51 * mrSges(7,2) - t100;
t101 = m(7) * pkin(5) + t51 * mrSges(7,1) - t47 * mrSges(7,2) + mrSges(6,1);
t71 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t87 = -mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t90 = m(5) + t96;
t82 = m(4) + t90;
t99 = t82 * qJ(3) - t87;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t98 = -t96 * qJ(3) - t101 * t52 + t71 * t48 + t87;
t45 = sin(pkin(6));
t50 = sin(qJ(1));
t95 = t45 * t50;
t54 = cos(qJ(1));
t94 = t45 * t54;
t53 = cos(qJ(2));
t93 = t50 * t53;
t92 = t54 * t53;
t91 = t54 * pkin(1) + pkin(8) * t95;
t89 = pkin(6) - qJ(2);
t88 = pkin(6) + qJ(2);
t72 = sin(t88);
t66 = t72 / 0.2e1;
t73 = sin(t89);
t62 = t66 - t73 / 0.2e1;
t27 = -t50 * t62 + t92;
t84 = t27 * pkin(2) + t91;
t80 = -t50 * pkin(1) + pkin(8) * t94;
t49 = sin(qJ(2));
t68 = cos(t89) / 0.2e1;
t74 = cos(t88);
t57 = t68 + t74 / 0.2e1;
t23 = t49 * t50 - t54 * t57;
t11 = t23 * pkin(2);
t67 = t73 / 0.2e1;
t63 = t67 - t72 / 0.2e1;
t25 = -t54 * t63 + t93;
t79 = qJ(3) * t25 - t11;
t26 = t54 * t49 + t50 * t57;
t17 = t26 * pkin(2);
t28 = t50 * t63 + t92;
t78 = qJ(3) * t28 - t17;
t35 = t66 + t67;
t33 = t35 * pkin(2);
t36 = t68 - t74 / 0.2e1;
t77 = qJ(3) * t36 + t33;
t75 = t27 * pkin(3) + t84;
t24 = t54 * t62 + t93;
t70 = -t24 * pkin(2) + t80;
t69 = t26 * pkin(4) + t75;
t64 = -t24 * pkin(3) + t70;
t3 = -t23 * t48 + t52 * t94;
t4 = t23 * t52 + t48 * t94;
t56 = mrSges(2,2) + (t90 * qJ(4) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3)) * t45;
t46 = cos(pkin(6));
t32 = t35 * pkin(3);
t22 = -t35 * t52 - t46 * t48;
t16 = t26 * pkin(3);
t10 = t23 * pkin(3);
t8 = t26 * t52 - t48 * t95;
t7 = t26 * t48 + t52 * t95;
t2 = t27 * t47 + t51 * t8;
t1 = t27 * t51 - t47 * t8;
t5 = [(-t54 * mrSges(2,1) - m(3) * t91 - m(4) * t84 - m(5) * t75 - m(6) * t69 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t69) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t71 * t7 - t99 * t26 + t56 * t50 + t100 * t27) * g(2) + (t50 * mrSges(2,1) - m(3) * t80 - m(4) * t70 - m(5) * t64 + t71 * t3 + t101 * t4 + t99 * t23 + t56 * t54 + t97 * t24 + t96 * (t23 * pkin(4) - t64)) * g(1) (-m(4) * t77 - m(5) * (t32 + t77) - t96 * (t36 * pkin(4) + t32 + t33) + t98 * t36 - t97 * t35) * g(3) + (-m(4) * t79 - m(5) * (-t10 + t79) - t96 * (t25 * pkin(4) - t10 - t11) + t98 * t25 + t97 * t23) * g(2) + (-m(4) * t78 - m(5) * (-t16 + t78) - t96 * (t28 * pkin(4) - t16 - t17) + t98 * t28 + t97 * t26) * g(1) (-g(1) * t26 - g(2) * t23 + g(3) * t35) * t82 (t46 * g(3) + (t50 * g(1) - t54 * g(2)) * t45) * t90 (t71 * t22 - t101 * (t35 * t48 - t46 * t52)) * g(3) + (-t101 * t3 + t71 * t4) * g(2) + (t101 * t7 + t71 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t24 * t51 - t4 * t47) * mrSges(7,1) + (-t24 * t47 - t4 * t51) * mrSges(7,2)) - g(3) * ((-t22 * t47 + t36 * t51) * mrSges(7,1) + (-t22 * t51 - t36 * t47) * mrSges(7,2))];
taug  = t5(:);
