% Calculate Gravitation load on the joints for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:04
% EndTime: 2018-11-23 16:24:05
% DurationCPUTime: 0.72s
% Computational Cost: add. (522->100), mult. (478->117), div. (0->0), fcn. (429->10), ass. (0->54)
t103 = -mrSges(6,1) - mrSges(7,1);
t101 = -mrSges(7,2) - mrSges(6,3);
t37 = sin(qJ(5));
t40 = cos(qJ(5));
t102 = -t37 * mrSges(7,3) + t103 * t40;
t97 = -m(6) - m(7);
t36 = qJ(3) + qJ(4);
t31 = sin(t36);
t65 = qJ(6) * t37;
t93 = (-m(7) * (-pkin(5) * t40 - pkin(4) - t65) - t102) * t31;
t32 = cos(t36);
t99 = -t32 * mrSges(5,1) + (mrSges(5,2) + t101) * t31;
t98 = m(3) + m(4);
t66 = t32 * pkin(4) + t31 * pkin(9);
t82 = pkin(4) * t31;
t38 = sin(qJ(3));
t83 = pkin(3) * t38;
t95 = m(7) * t83 - m(6) * (-t82 - t83) + t93;
t35 = qJ(1) + pkin(10);
t30 = cos(t35);
t74 = t30 * t32;
t14 = pkin(9) * t74;
t72 = t31 * t37;
t64 = mrSges(6,2) * t72;
t94 = -m(7) * t14 + t101 * t74 - t30 * t64;
t29 = sin(t35);
t92 = g(1) * t30 + g(2) * t29;
t91 = -m(5) + t97;
t71 = t32 * t37;
t90 = mrSges(6,2) * t71 + t102 * t32 + t99;
t89 = m(7) * pkin(5) - t103;
t88 = (-t64 + (pkin(9) * t97 + t101) * t32) * t29;
t87 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t41 = cos(qJ(3));
t55 = t41 * mrSges(4,1) - t38 * mrSges(4,2);
t86 = -m(4) * pkin(2) - mrSges(3,1) - t55 + t99;
t85 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t39 = sin(qJ(1));
t78 = t39 * pkin(1);
t33 = t41 * pkin(3);
t42 = cos(qJ(1));
t34 = t42 * pkin(1);
t75 = t30 * t31;
t70 = t32 * t40;
t58 = pkin(5) * t70 + t32 * t65 + t66;
t28 = t33 + pkin(2);
t43 = -pkin(8) - pkin(7);
t57 = t30 * t28 - t29 * t43 + t34;
t52 = mrSges(5,1) * t31 + mrSges(5,2) * t32;
t4 = t29 * t37 + t30 * t70;
t3 = -t29 * t40 + t30 * t71;
t2 = t29 * t70 - t30 * t37;
t1 = t29 * t71 + t30 * t40;
t5 = [(-m(5) * t57 - t42 * mrSges(2,1) + t39 * mrSges(2,2) + t97 * (pkin(4) * t74 + pkin(9) * t75 + t57) - t89 * t4 - t98 * t34 - t87 * t3 + t85 * t29 + t86 * t30) * g(2) + (t39 * mrSges(2,1) + t42 * mrSges(2,2) + t98 * t78 + t91 * (-t30 * t43 - t78) + t89 * t2 + t87 * t1 + t85 * t30 + (m(5) * t28 + t97 * (-t28 - t66) - t86) * t29) * g(1) (t91 - t98) * g(3) (t95 * t29 + t88) * g(2) + (-m(6) * t14 + t95 * t30 + t94) * g(1) + (-t55 - m(5) * t33 - m(6) * (t33 + t66) - m(7) * (t33 + t58) + t90) * g(3) + (m(5) * t83 + mrSges(4,1) * t38 + mrSges(4,2) * t41 + t52) * t92, t92 * t52 + ((m(6) * t82 + t93) * t29 + t88) * g(2) + (-m(6) * (-pkin(4) * t75 + t14) + t93 * t30 + t94) * g(1) + (-m(6) * t66 - m(7) * t58 + t90) * g(3) (t89 * t37 - t87 * t40) * g(3) * t31 + (t89 * t1 - t87 * t2) * g(2) + (t89 * t3 - t87 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t72) * m(7)];
taug  = t5(:);
