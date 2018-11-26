% Calculate Gravitation load on the joints for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:02:35
% EndTime: 2018-11-23 15:02:36
% DurationCPUTime: 0.70s
% Computational Cost: add. (1015->96), mult. (1223->136), div. (0->0), fcn. (1185->14), ass. (0->61)
t100 = -m(4) - m(5);
t99 = -m(6) - m(7);
t96 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t73 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t69 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t98 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t97 = t69 * t50 - t73 * t53 - mrSges(5,1);
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t95 = -t51 * mrSges(5,1) + mrSges(3,2) - mrSges(4,3) + t99 * (-pkin(9) * t54 + qJ(3)) + t100 * qJ(3) - t96 * t54;
t94 = pkin(4) * t51;
t47 = sin(pkin(10));
t48 = sin(pkin(6));
t93 = t47 * t48;
t55 = cos(qJ(2));
t92 = t47 * t55;
t91 = t50 * t51;
t90 = t51 * t53;
t89 = cos(pkin(10));
t88 = pkin(6) - qJ(2);
t87 = pkin(6) + qJ(2);
t86 = -t99 - t100;
t52 = sin(qJ(2));
t68 = cos(t88) / 0.2e1;
t76 = cos(t87);
t64 = t68 + t76 / 0.2e1;
t29 = t47 * t52 - t89 * t64;
t26 = t29 * pkin(2);
t85 = -pkin(8) * t29 - t26;
t31 = t47 * t64 + t89 * t52;
t27 = t31 * pkin(2);
t84 = -pkin(8) * t31 - t27;
t74 = sin(t87);
t66 = t74 / 0.2e1;
t75 = sin(t88);
t67 = t75 / 0.2e1;
t41 = t66 + t67;
t40 = t41 * pkin(2);
t83 = pkin(8) * t41 + t40;
t78 = t48 * t89;
t77 = t89 * t55;
t63 = t67 - t74 / 0.2e1;
t62 = t66 - t75 / 0.2e1;
t57 = t89 * t62 + t92;
t56 = -t47 * t62 + t77;
t49 = cos(pkin(6));
t42 = t68 - t76 / 0.2e1;
t34 = -t41 * t51 + t49 * t54;
t33 = -t41 * t54 - t49 * t51;
t32 = t47 * t63 + t77;
t30 = -t89 * t63 + t92;
t16 = -t29 * t51 + t54 * t78;
t15 = t29 * t54 + t51 * t78;
t14 = t31 * t51 + t54 * t93;
t13 = t31 * t54 - t51 * t93;
t9 = t34 * t50 - t42 * t53;
t3 = -t16 * t50 - t57 * t53;
t1 = t14 * t50 - t56 * t53;
t2 = [(-m(2) - m(3) - t86) * g(3) (-m(4) * t40 - m(5) * t83 + t99 * (t42 * t94 + t83) - t73 * (t41 * t50 + t42 * t90) + t69 * (-t41 * t53 + t42 * t91) - t98 * t41 + t95 * t42) * g(3) + (m(4) * t26 - m(5) * t85 + t99 * (t30 * t94 + t85) - t73 * (-t29 * t50 + t30 * t90) + t69 * (t29 * t53 + t30 * t91) + t98 * t29 + t95 * t30) * g(2) + (m(4) * t27 - m(5) * t84 + t99 * (t32 * t94 + t84) - t73 * (-t31 * t50 + t32 * t90) + t69 * (t31 * t53 + t32 * t91) + t98 * t31 + t95 * t32) * g(1) (-g(1) * t31 - g(2) * t29 + g(3) * t41) * t86 (t99 * (t33 * pkin(4) + pkin(9) * t34) + t96 * t34 + t97 * t33) * g(3) + (t99 * (t15 * pkin(4) - pkin(9) * t16) - t96 * t16 + t97 * t15) * g(2) + (t99 * (t13 * pkin(4) + pkin(9) * t14) + t96 * t14 + t97 * t13) * g(1) (t73 * t9 + t69 * (t34 * t53 + t42 * t50)) * g(3) + (t69 * (-t16 * t53 + t57 * t50) + t73 * t3) * g(2) + (t69 * (t14 * t53 + t56 * t50) + t73 * t1) * g(1) (-g(1) * t1 - g(2) * t3 - g(3) * t9) * m(7)];
taug  = t2(:);
