% Calculate Gravitation load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:45
% EndTime: 2018-11-23 15:13:46
% DurationCPUTime: 0.84s
% Computational Cost: add. (1168->106), mult. (1433->149), div. (0->0), fcn. (1410->14), ass. (0->64)
t115 = -m(6) - m(7);
t114 = m(5) - t115;
t117 = -mrSges(4,2) + mrSges(5,3);
t112 = mrSges(4,1) - mrSges(5,2) + mrSges(7,2) + mrSges(6,3);
t84 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t116 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t113 = mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t111 = t112 * t63 + t117 * t60 + mrSges(3,1);
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t110 = -t114 * qJ(4) + t116 * t62 - t84 * t59 - t117;
t58 = sin(pkin(10));
t61 = sin(qJ(2));
t96 = pkin(6) + qJ(2);
t82 = cos(t96) / 0.2e1;
t97 = pkin(6) - qJ(2);
t87 = cos(t97);
t67 = t87 / 0.2e1 + t82;
t99 = cos(pkin(10));
t35 = t58 * t61 - t67 * t99;
t109 = t35 * t63;
t38 = t58 * t67 + t61 * t99;
t108 = t38 * t63;
t85 = sin(t96);
t80 = t85 / 0.2e1;
t86 = sin(t97);
t81 = t86 / 0.2e1;
t51 = t80 + t81;
t107 = t51 * t63;
t64 = cos(qJ(2));
t106 = t58 * t64;
t105 = t59 * t60;
t104 = t60 * t62;
t101 = qJ(4) * t60;
t100 = cos(pkin(6));
t98 = sin(pkin(6));
t66 = t81 - t85 / 0.2e1;
t37 = -t66 * t99 + t106;
t95 = -t35 * pkin(2) + pkin(8) * t37;
t88 = t99 * t64;
t40 = t58 * t66 + t88;
t94 = -t38 * pkin(2) + pkin(8) * t40;
t52 = t82 - t87 / 0.2e1;
t93 = t51 * pkin(2) - pkin(8) * t52;
t89 = t58 * t98;
t77 = -pkin(3) * t109 - t35 * t101 + t95;
t76 = -pkin(3) * t108 - t38 * t101 + t94;
t75 = pkin(3) * t107 + t51 * t101 + t93;
t74 = t99 * t98;
t73 = t80 - t86 / 0.2e1;
t41 = -t100 * t63 - t52 * t60;
t39 = -t58 * t73 + t88;
t36 = t73 * t99 + t106;
t34 = t41 * pkin(3);
t17 = t39 * t60 - t63 * t89;
t15 = t36 * t60 + t63 * t74;
t14 = t17 * pkin(3);
t13 = t15 * pkin(3);
t9 = -t41 * t62 - t51 * t59;
t3 = -t17 * t62 + t38 * t59;
t1 = -t15 * t62 + t35 * t59;
t2 = [(-m(2) - m(3) - m(4) - t114) * g(3) (-m(4) * t93 - m(5) * t75 + t115 * (-t52 * pkin(4) + pkin(9) * t107 + t75) - t84 * (t105 * t51 - t52 * t62) - t116 * (-t104 * t51 - t52 * t59) - t113 * t52 - t111 * t51) * g(3) + (-m(4) * t95 - m(5) * t77 + t115 * (t37 * pkin(4) - pkin(9) * t109 + t77) - t84 * (-t105 * t35 + t37 * t62) - t116 * (t104 * t35 + t37 * t59) + t113 * t37 + t111 * t35) * g(2) + (-m(4) * t94 - m(5) * t76 + t115 * (t40 * pkin(4) - pkin(9) * t108 + t76) - t84 * (-t105 * t38 + t40 * t62) - t116 * (t104 * t38 + t40 * t59) + t113 * t40 + t111 * t38) * g(1) (m(5) * t34 + t115 * (-pkin(9) * t41 - t34) + t110 * (t100 * t60 - t52 * t63) + t112 * t41) * g(3) + (m(5) * t13 + t115 * (-pkin(9) * t15 - t13) + t110 * (t36 * t63 - t60 * t74) + t112 * t15) * g(2) + (m(5) * t14 + t115 * (-pkin(9) * t17 - t14) + t110 * (t39 * t63 + t60 * t89) + t112 * t17) * g(1), t114 * (-g(1) * t17 - g(2) * t15 - g(3) * t41) (t84 * t9 - t116 * (t41 * t59 - t51 * t62)) * g(3) + (-t116 * (t15 * t59 + t35 * t62) + t84 * t1) * g(2) + (-t116 * (t17 * t59 + t38 * t62) + t84 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
