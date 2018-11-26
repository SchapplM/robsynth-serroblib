% Calculate Gravitation load on the joints for
% S6PPRRRP2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:51:23
% EndTime: 2018-11-23 14:51:23
% DurationCPUTime: 0.80s
% Computational Cost: add. (3165->96), mult. (3319->141), div. (0->0), fcn. (3280->22), ass. (0->79)
t102 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t99 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t133 = -t102 * t72 + t99 * t69 - mrSges(5,1);
t132 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t136 = -m(6) - m(7);
t142 = -m(5) + t136;
t115 = pkin(7) - qJ(3);
t103 = sin(t115);
t114 = pkin(7) + qJ(3);
t96 = sin(t114) / 0.2e1;
t134 = t96 - t103 / 0.2e1;
t112 = pkin(6) + pkin(12);
t101 = cos(t112);
t113 = pkin(6) - pkin(12);
t90 = cos(t113) / 0.2e1;
t62 = t90 - t101 / 0.2e1;
t74 = cos(qJ(3));
t100 = sin(t113);
t89 = sin(t112) / 0.2e1;
t80 = t89 + t100 / 0.2e1;
t141 = t134 * t80 + t62 * t74;
t117 = sin(pkin(11));
t118 = cos(pkin(11));
t61 = t89 - t100 / 0.2e1;
t67 = cos(pkin(12));
t55 = -t117 * t61 + t118 * t67;
t116 = sin(pkin(12));
t81 = t90 + t101 / 0.2e1;
t76 = t116 * t118 + t117 * t81;
t140 = -t134 * t76 + t55 * t74;
t54 = t117 * t67 + t118 * t61;
t75 = t116 * t117 - t118 * t81;
t139 = -t134 * t75 + t54 * t74;
t138 = -t102 * t69 - t99 * t72 + mrSges(4,2) - mrSges(5,3);
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t137 = mrSges(4,1) + (-t136 * pkin(10) - t132) * t70 + (-t136 * pkin(4) - t133) * t73;
t131 = m(3) + m(4) - t142;
t66 = sin(pkin(6));
t68 = cos(pkin(7));
t124 = t66 * t68;
t119 = cos(pkin(6));
t105 = cos(t115);
t104 = cos(t114);
t98 = t105 / 0.2e1;
t97 = t104 / 0.2e1;
t88 = t98 - t104 / 0.2e1;
t86 = t66 * t88;
t85 = t98 + t97;
t84 = t97 - t105 / 0.2e1;
t82 = t96 + t103 / 0.2e1;
t79 = t66 * t84;
t78 = t66 * t82;
t71 = sin(qJ(3));
t65 = sin(pkin(7));
t56 = t119 * t68 - t65 * t80;
t47 = t117 * t124 + t65 * t76;
t46 = -t118 * t124 + t65 * t75;
t45 = -t119 * t84 + t141;
t44 = t119 * t88 + t141;
t43 = -t119 * t82 + t62 * t71 - t80 * t85;
t37 = -t117 * t79 + t140;
t36 = t117 * t86 + t140;
t35 = -t117 * t78 + t55 * t71 + t76 * t85;
t34 = t118 * t79 + t139;
t33 = -t118 * t86 + t139;
t32 = t118 * t78 + t54 * t71 + t75 * t85;
t21 = t44 * t73 + t56 * t70;
t20 = -t44 * t70 + t56 * t73;
t16 = t36 * t73 + t47 * t70;
t15 = -t36 * t70 + t47 * t73;
t14 = t33 * t73 + t46 * t70;
t13 = -t33 * t70 + t46 * t73;
t9 = t21 * t69 - t43 * t72;
t3 = t16 * t69 - t35 * t72;
t1 = t14 * t69 - t32 * t72;
t2 = [(-m(2) - t131) * g(3) (-t119 * g(3) + (-g(1) * t117 + g(2) * t118) * t66) * t131 (t142 * (-t43 * pkin(3) + t45 * pkin(9)) + t138 * t45 + t137 * t43) * g(3) + (t142 * (-t32 * pkin(3) + pkin(9) * t34) + t138 * t34 + t137 * t32) * g(2) + (t142 * (-t35 * pkin(3) + pkin(9) * t37) + t138 * t37 + t137 * t35) * g(1) (t136 * (t20 * pkin(4) + pkin(10) * t21) + t132 * t21 + t133 * t20) * g(3) + (t136 * (t13 * pkin(4) + pkin(10) * t14) + t132 * t14 + t133 * t13) * g(2) + (t136 * (t15 * pkin(4) + pkin(10) * t16) + t132 * t16 + t133 * t15) * g(1) (t102 * t9 + t99 * (t21 * t72 + t43 * t69)) * g(3) + (t99 * (t14 * t72 + t32 * t69) + t102 * t1) * g(2) + (t99 * (t16 * t72 + t35 * t69) + t102 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
