% Calculate Gravitation load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:56:06
% EndTime: 2018-11-23 16:56:07
% DurationCPUTime: 0.85s
% Computational Cost: add. (1162->118), mult. (1249->148), div. (0->0), fcn. (1180->16), ass. (0->67)
t108 = m(4) + m(5);
t107 = m(6) + m(7);
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t106 = m(7) * pkin(5) + t60 * mrSges(7,1) - t57 * mrSges(7,2) + mrSges(6,1);
t110 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t109 = -t52 * mrSges(5,1) - t54 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t72 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t101 = mrSges(7,1) * t57 + mrSges(7,2) * t60 - t72;
t104 = m(5) + t107;
t88 = m(4) + t104;
t103 = t88 * qJ(3) - t109;
t51 = pkin(11) + qJ(5);
t48 = sin(t51);
t49 = cos(t51);
t102 = t109 - t106 * t48 + t110 * t49 + (-t107 - t108) * qJ(3);
t100 = pkin(4) * t52;
t53 = sin(pkin(6));
t59 = sin(qJ(1));
t99 = t53 * t59;
t62 = cos(qJ(1));
t98 = t53 * t62;
t61 = cos(qJ(2));
t97 = t59 * t61;
t96 = t62 * t61;
t95 = t62 * pkin(1) + pkin(8) * t99;
t94 = pkin(6) - qJ(2);
t93 = pkin(6) + qJ(2);
t83 = sin(t93);
t78 = t83 / 0.2e1;
t84 = sin(t94);
t73 = t78 - t84 / 0.2e1;
t30 = -t59 * t73 + t96;
t90 = t30 * pkin(2) + t95;
t86 = -t59 * pkin(1) + pkin(8) * t98;
t85 = cos(t93);
t27 = t62 * t73 + t97;
t81 = t27 * pkin(2) - t86;
t80 = cos(t94) / 0.2e1;
t79 = t84 / 0.2e1;
t58 = sin(qJ(2));
t65 = t80 + t85 / 0.2e1;
t29 = t62 * t58 + t59 * t65;
t47 = pkin(4) * t54 + pkin(3);
t56 = -pkin(9) - qJ(4);
t75 = t29 * t100 - t30 * t56 + t47 * t99 + t90;
t74 = t79 - t83 / 0.2e1;
t26 = t58 * t59 - t62 * t65;
t7 = -t26 * t48 + t49 * t98;
t5 = t26 * t49 + t48 * t98;
t64 = mrSges(2,2) + (-m(5) * pkin(3) - mrSges(5,1) * t54 + mrSges(5,2) * t52 - mrSges(4,1) - mrSges(3,3)) * t53;
t55 = cos(pkin(6));
t41 = t80 - t85 / 0.2e1;
t40 = t78 + t79;
t37 = t40 * pkin(2);
t31 = t59 * t74 + t96;
t28 = -t62 * t74 + t97;
t24 = t29 * pkin(2);
t22 = t26 * pkin(2);
t13 = -t40 * t48 + t49 * t55;
t4 = t29 * t48 + t49 * t99;
t3 = -t29 * t49 + t48 * t99;
t2 = t30 * t57 + t4 * t60;
t1 = t30 * t60 - t4 * t57;
t6 = [(-t62 * mrSges(2,1) - m(3) * t95 - m(6) * t75 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t75) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t110 * t3 - t103 * t29 + t64 * t59 + t72 * t30 - t108 * t90) * g(2) + (t59 * mrSges(2,1) - m(3) * t86 - t110 * t5 - t106 * t7 + t103 * t26 + t64 * t62 + t101 * t27 + t108 * t81 + t107 * (t26 * t100 - t27 * t56 - t47 * t98 + t81)) * g(1) (-t107 * (t41 * t100 - t40 * t56 + t37) - t108 * t37 + t102 * t41 - t101 * t40) * g(3) + (-t107 * (t28 * t100 + t26 * t56 - t22) + t108 * t22 + t102 * t28 + t101 * t26) * g(2) + (-t107 * (t31 * t100 + t29 * t56 - t24) + t108 * t24 + t102 * t31 + t101 * t29) * g(1) (-g(1) * t29 - g(2) * t26 + g(3) * t40) * t88, t104 * (-g(1) * t30 - g(2) * t27 - g(3) * t41) (-t110 * t13 - t106 * (-t40 * t49 - t48 * t55)) * g(3) + (-t106 * t5 + t110 * t7) * g(2) + (t106 * t3 - t110 * t4) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t27 * t60 + t57 * t7) * mrSges(7,1) + (-t27 * t57 + t60 * t7) * mrSges(7,2)) - g(3) * ((-t13 * t57 + t41 * t60) * mrSges(7,1) + (-t13 * t60 - t41 * t57) * mrSges(7,2))];
taug  = t6(:);
