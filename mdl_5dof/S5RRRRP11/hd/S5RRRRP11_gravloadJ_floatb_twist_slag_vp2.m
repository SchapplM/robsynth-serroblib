% Calculate Gravitation load on the joints for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:26
% EndTime: 2019-12-31 22:14:29
% DurationCPUTime: 0.89s
% Computational Cost: add. (509->103), mult. (1263->148), div. (0->0), fcn. (1511->10), ass. (0->55)
t118 = m(5) + m(6);
t126 = t118 * pkin(9);
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t79 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t80 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t125 = t79 * t60 - t80 * t63;
t124 = m(4) + t118;
t123 = mrSges(5,3) + mrSges(6,2);
t107 = sin(qJ(1));
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t108 = cos(qJ(1));
t94 = cos(pkin(5));
t77 = t94 * t108;
t42 = t107 * t65 + t62 * t77;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t59 = sin(pkin(5));
t89 = t59 * t108;
t18 = t42 * t64 - t61 * t89;
t41 = t107 * t62 - t65 * t77;
t1 = t18 * t60 - t41 * t63;
t121 = t18 * t63 + t41 * t60;
t97 = mrSges(4,3) - mrSges(3,2);
t120 = -t80 * t60 - t79 * t63 - t97;
t117 = -t64 * mrSges(4,1) + t61 * mrSges(4,2) - mrSges(3,1);
t119 = -t117 + (t123 + t126) * t61 + (t118 * pkin(3) - t125) * t64;
t115 = -mrSges(4,1) + t125;
t114 = mrSges(4,2) - t123;
t113 = pkin(8) * t124 + t97;
t104 = t59 * t62;
t103 = t59 * t65;
t98 = t64 * t65;
t96 = pkin(2) * t103 + pkin(8) * t104;
t88 = t59 * t107;
t95 = t108 * pkin(1) + pkin(7) * t88;
t92 = t61 * t103;
t91 = t59 * t98;
t76 = t94 * t107;
t44 = t108 * t65 - t62 * t76;
t90 = t44 * pkin(2) + t95;
t17 = -t42 * t61 - t64 * t89;
t78 = -pkin(1) * t107 + pkin(7) * t89;
t69 = -t42 * pkin(2) + t78;
t68 = t114 - t126;
t43 = t108 * t62 + t65 * t76;
t40 = t104 * t64 + t61 * t94;
t39 = -t104 * t61 + t64 * t94;
t22 = t44 * t64 + t61 * t88;
t21 = t44 * t61 - t64 * t88;
t15 = t103 * t63 + t40 * t60;
t6 = t22 * t63 + t43 * t60;
t5 = t22 * t60 - t43 * t63;
t2 = [(-t108 * mrSges(2,1) + t107 * mrSges(2,2) - m(3) * t95 - t44 * mrSges(3,1) - mrSges(3,3) * t88 - m(4) * t90 - t22 * mrSges(4,1) - t80 * t6 + t79 * t5 - t113 * t43 + t68 * t21 - t118 * (t22 * pkin(3) + t90)) * g(2) + (t107 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t78 + t42 * mrSges(3,1) - mrSges(3,3) * t89 - m(4) * t69 + t18 * mrSges(4,1) + t113 * t41 + t80 * t121 - t79 * t1 + t68 * t17 + t118 * (pkin(3) * t18 - t69)) * g(1), (-m(4) * t96 - t123 * t92 - t118 * (pkin(3) * t91 + pkin(9) * t92 + t96) + t79 * (-t104 * t63 + t60 * t91) + (t117 * t65 - t97 * t62 - t80 * (t60 * t62 + t63 * t98)) * t59) * g(3) + (-t124 * (-t41 * pkin(2) + pkin(8) * t42) + t120 * t42 + t119 * t41) * g(2) + (-t124 * (-t43 * pkin(2) + pkin(8) * t44) + t120 * t44 + t119 * t43) * g(1), (-t118 * (t39 * pkin(3) + pkin(9) * t40) + t114 * t40 + t115 * t39) * g(3) + (-t118 * (t17 * pkin(3) + pkin(9) * t18) + t114 * t18 + t115 * t17) * g(2) + (-t118 * (-t21 * pkin(3) + pkin(9) * t22) + t114 * t22 - t115 * t21) * g(1), (t79 * (-t103 * t60 + t40 * t63) + t80 * t15) * g(3) + (t80 * t1 + t79 * t121) * g(2) + (t5 * t80 + t6 * t79) * g(1), (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(6)];
taug = t2(:);
