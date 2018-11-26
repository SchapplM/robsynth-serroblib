% Calculate Gravitation load on the joints for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2018-11-23 15:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:01:24
% EndTime: 2018-11-23 15:01:25
% DurationCPUTime: 0.68s
% Computational Cost: add. (1173->100), mult. (1212->141), div. (0->0), fcn. (1176->16), ass. (0->58)
t79 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t78 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t109 = -m(6) - m(7);
t108 = mrSges(6,3) + mrSges(7,2);
t58 = pkin(11) + qJ(4);
t56 = sin(t58);
t57 = cos(t58);
t62 = cos(pkin(11));
t107 = mrSges(3,1) + m(4) * pkin(2) + t62 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + t57 * mrSges(5,1) - mrSges(5,2) * t56;
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t106 = t78 * t65 - t79 * t67 - mrSges(5,1);
t105 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t104 = mrSges(5,2) - t108;
t103 = pkin(4) * t57;
t60 = sin(pkin(10));
t66 = sin(qJ(2));
t91 = pkin(6) + qJ(2);
t77 = cos(t91) / 0.2e1;
t92 = pkin(6) - qJ(2);
t81 = cos(t92);
t70 = t81 / 0.2e1 + t77;
t93 = cos(pkin(10));
t34 = t60 * t66 - t93 * t70;
t102 = t34 * t56;
t37 = t60 * t70 + t93 * t66;
t101 = t37 * t56;
t76 = sin(t91) / 0.2e1;
t80 = sin(t92);
t48 = t76 + t80 / 0.2e1;
t100 = t48 * t56;
t99 = t57 * t65;
t98 = t57 * t67;
t61 = sin(pkin(6));
t97 = t60 * t61;
t55 = pkin(3) * t62 + pkin(2);
t64 = -pkin(8) - qJ(3);
t49 = t76 - t80 / 0.2e1;
t68 = cos(qJ(2));
t73 = t93 * t49 + t60 * t68;
t96 = -t34 * t55 - t64 * t73;
t72 = -t60 * t49 + t93 * t68;
t95 = -t37 * t55 - t64 * t72;
t50 = t77 - t81 / 0.2e1;
t94 = t48 * t55 + t50 * t64;
t90 = m(4) + m(5) - t109;
t85 = t61 * t93;
t63 = cos(pkin(6));
t27 = -t50 * t57 + t56 * t63;
t26 = t50 * t56 + t57 * t63;
t16 = t56 * t97 + t57 * t72;
t15 = -t56 * t72 + t57 * t97;
t14 = -t56 * t85 + t57 * t73;
t13 = -t56 * t73 - t57 * t85;
t9 = t27 * t65 + t48 * t67;
t3 = t16 * t65 - t37 * t67;
t1 = t14 * t65 - t34 * t67;
t2 = [(-m(2) - m(3) - t90) * g(3) (-m(5) * t94 + t109 * (pkin(9) * t100 + t48 * t103 + t94) - t79 * (t48 * t98 - t50 * t65) + t78 * (t48 * t99 + t50 * t67) - t108 * t100 - t105 * t50 - t107 * t48) * g(3) + (-m(5) * t96 + t109 * (-pkin(9) * t102 - t34 * t103 + t96) - t79 * (-t34 * t98 + t65 * t73) + t78 * (-t34 * t99 - t67 * t73) + t108 * t102 + t105 * t73 + t107 * t34) * g(2) + (-m(5) * t95 + t109 * (-pkin(9) * t101 - t37 * t103 + t95) - t79 * (-t37 * t98 + t65 * t72) + t78 * (-t37 * t99 - t67 * t72) + t108 * t101 + t105 * t72 + t107 * t37) * g(1) (-g(1) * t37 - g(2) * t34 + g(3) * t48) * t90 (t109 * (t26 * pkin(4) + pkin(9) * t27) + t104 * t27 + t106 * t26) * g(3) + (t109 * (t13 * pkin(4) + pkin(9) * t14) + t104 * t14 + t106 * t13) * g(2) + (t109 * (t15 * pkin(4) + pkin(9) * t16) + t104 * t16 + t106 * t15) * g(1) (t79 * t9 + t78 * (t27 * t67 - t48 * t65)) * g(3) + (t78 * (t14 * t67 + t34 * t65) + t79 * t1) * g(2) + (t78 * (t16 * t67 + t37 * t65) + t79 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
