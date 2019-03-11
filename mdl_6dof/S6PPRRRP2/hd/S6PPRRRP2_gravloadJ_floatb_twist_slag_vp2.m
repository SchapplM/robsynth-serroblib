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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:03
% EndTime: 2019-03-08 18:56:05
% DurationCPUTime: 0.70s
% Computational Cost: add. (924->97), mult. (2572->157), div. (0->0), fcn. (3280->14), ass. (0->66)
t82 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t81 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t113 = -m(6) - m(7);
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t112 = t58 * mrSges(5,1) - mrSges(5,2) * t55 + mrSges(4,1);
t111 = mrSges(4,2) - mrSges(5,3);
t110 = mrSges(6,3) + mrSges(7,2);
t89 = sin(pkin(12));
t90 = sin(pkin(11));
t67 = t90 * t89;
t93 = cos(pkin(12));
t94 = cos(pkin(11));
t74 = t94 * t93;
t96 = cos(pkin(6));
t61 = -t96 * t74 + t67;
t91 = sin(pkin(7));
t92 = sin(pkin(6));
t71 = t92 * t91;
t95 = cos(pkin(7));
t109 = t61 * t95 + t94 * t71;
t68 = t90 * t93;
t72 = t94 * t89;
t62 = t96 * t68 + t72;
t70 = t92 * t90;
t108 = t62 * t95 - t91 * t70;
t107 = t93 * t95 * t92 + t96 * t91;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t106 = t81 * t54 - t82 * t57 - mrSges(5,1);
t105 = mrSges(5,2) - t110;
t104 = m(3) + m(4) + m(5) - t113;
t103 = pkin(4) * t58;
t102 = cos(qJ(3));
t48 = t96 * t72 + t68;
t56 = sin(qJ(3));
t30 = t109 * t102 + t48 * t56;
t101 = t30 * t55;
t49 = -t96 * t67 + t74;
t32 = t108 * t102 + t49 * t56;
t100 = t32 * t55;
t69 = t92 * t89;
t41 = -t107 * t102 + t56 * t69;
t99 = t41 * t55;
t98 = t54 * t58;
t97 = t57 * t58;
t31 = t48 * t102 - t109 * t56;
t87 = -t30 * pkin(3) + pkin(9) * t31;
t33 = t49 * t102 - t108 * t56;
t86 = -t32 * pkin(3) + pkin(9) * t33;
t42 = t102 * t69 + t107 * t56;
t85 = -t41 * pkin(3) + pkin(9) * t42;
t73 = t94 * t92;
t47 = -t93 * t71 + t96 * t95;
t44 = t62 * t91 + t95 * t70;
t43 = t61 * t91 - t95 * t73;
t35 = t42 * t58 + t47 * t55;
t34 = -t42 * t55 + t47 * t58;
t15 = t35 * t54 - t41 * t57;
t14 = t33 * t58 + t44 * t55;
t13 = -t33 * t55 + t44 * t58;
t12 = t31 * t58 + t43 * t55;
t11 = -t31 * t55 + t43 * t58;
t3 = t14 * t54 - t32 * t57;
t1 = t12 * t54 - t30 * t57;
t2 = [(-m(2) - t104) * g(3) (-g(1) * t70 + g(2) * t73 - g(3) * t96) * t104 (-m(5) * t85 + t110 * t99 + t113 * (-pkin(10) * t99 - t41 * t103 + t85) + t111 * t42 + t112 * t41 - t82 * (-t41 * t97 + t42 * t54) + t81 * (-t41 * t98 - t42 * t57)) * g(3) + (-m(5) * t87 + t113 * (-pkin(10) * t101 - t30 * t103 + t87) - t82 * (-t30 * t97 + t31 * t54) + t81 * (-t30 * t98 - t31 * t57) + t111 * t31 + t112 * t30 + t110 * t101) * g(2) + (-m(5) * t86 + t113 * (-pkin(10) * t100 - t32 * t103 + t86) - t82 * (-t32 * t97 + t33 * t54) + t81 * (-t32 * t98 - t33 * t57) + t111 * t33 + t112 * t32 + t110 * t100) * g(1) (t113 * (t34 * pkin(4) + pkin(10) * t35) + t105 * t35 + t106 * t34) * g(3) + (t113 * (t11 * pkin(4) + pkin(10) * t12) + t105 * t12 + t106 * t11) * g(2) + (t113 * (t13 * pkin(4) + pkin(10) * t14) + t105 * t14 + t106 * t13) * g(1) (t81 * (t35 * t57 + t41 * t54) + t82 * t15) * g(3) + (t81 * (t12 * t57 + t30 * t54) + t82 * t1) * g(2) + (t81 * (t14 * t57 + t32 * t54) + t82 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
