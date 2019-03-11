% Calculate Gravitation load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:41
% EndTime: 2019-03-08 23:15:43
% DurationCPUTime: 1.18s
% Computational Cost: add. (636->108), mult. (1216->151), div. (0->0), fcn. (1417->14), ass. (0->50)
t38 = qJ(4) + pkin(12);
t34 = cos(t38);
t44 = cos(qJ(4));
t36 = t44 * pkin(4);
t21 = pkin(5) * t34 + t36;
t35 = qJ(6) + t38;
t30 = sin(t35);
t31 = cos(t35);
t33 = sin(t38);
t41 = sin(qJ(4));
t86 = mrSges(4,1) + m(7) * (pkin(3) + t21) + t31 * mrSges(7,1) - t30 * mrSges(7,2) + m(6) * (t36 + pkin(3)) + t34 * mrSges(6,1) - t33 * mrSges(6,2) + m(5) * pkin(3) + t44 * mrSges(5,1) - t41 * mrSges(5,2);
t40 = -qJ(5) - pkin(9);
t85 = mrSges(4,2) + m(7) * (-pkin(10) + t40) - mrSges(7,3) + m(6) * t40 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t101 = t85 * t42 - t86 * t45 - mrSges(3,1);
t93 = m(6) + m(7);
t87 = m(4) + m(5) + t93;
t77 = pkin(4) * t41;
t20 = pkin(5) * t33 + t77;
t96 = -m(6) * t77 - m(7) * t20 - t41 * mrSges(5,1) - t33 * mrSges(6,1) - t30 * mrSges(7,1) - t44 * mrSges(5,2) - t34 * mrSges(6,2) - t31 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t95 = pkin(2) * t87 - t101;
t91 = -m(6) * pkin(4) - mrSges(5,1);
t82 = -pkin(8) * t87 + t96;
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t72 = cos(pkin(11));
t73 = cos(pkin(6));
t59 = t73 * t72;
t71 = sin(pkin(11));
t13 = t43 * t71 - t46 * t59;
t14 = t43 * t59 + t46 * t71;
t39 = sin(pkin(6));
t67 = t39 * t72;
t8 = t14 * t45 - t42 * t67;
t80 = (t13 * t31 - t30 * t8) * mrSges(7,1) + (-t13 * t30 - t31 * t8) * mrSges(7,2);
t58 = t73 * t71;
t16 = -t43 * t58 + t46 * t72;
t66 = t39 * t71;
t10 = t16 * t45 + t42 * t66;
t15 = t43 * t72 + t46 * t58;
t79 = (-t10 * t30 + t15 * t31) * mrSges(7,1) + (-t10 * t31 - t15 * t30) * mrSges(7,2);
t76 = t39 * t43;
t18 = t42 * t73 + t45 * t76;
t75 = t39 * t46;
t78 = (-t18 * t30 - t31 * t75) * mrSges(7,1) + (-t18 * t31 + t30 * t75) * mrSges(7,2);
t17 = t42 * t76 - t45 * t73;
t9 = t16 * t42 - t45 * t66;
t7 = t14 * t42 + t45 * t67;
t1 = [(-m(2) - m(3) - t87) * g(3) (-t87 * (pkin(2) * t75 + pkin(8) * t76) + (t101 * t46 + t96 * t43) * t39) * g(3) + (t95 * t13 + t82 * t14) * g(2) + (t95 * t15 + t82 * t16) * g(1) (t86 * t17 + t85 * t18) * g(3) + (t86 * t7 + t85 * t8) * g(2) + (t85 * t10 + t86 * t9) * g(1) (-(-t18 * t44 + t41 * t75) * mrSges(5,2) - (-t18 * t33 - t34 * t75) * mrSges(6,1) - (-t18 * t34 + t33 * t75) * mrSges(6,2) - m(7) * (-t18 * t20 - t21 * t75) - t78 + t91 * (-t18 * t41 - t44 * t75)) * g(3) + (-(-t13 * t41 - t44 * t8) * mrSges(5,2) - (t13 * t34 - t33 * t8) * mrSges(6,1) - (-t13 * t33 - t34 * t8) * mrSges(6,2) - m(7) * (t13 * t21 - t20 * t8) - t80 + t91 * (t13 * t44 - t41 * t8)) * g(2) + (-(-t10 * t44 - t15 * t41) * mrSges(5,2) - (-t10 * t33 + t15 * t34) * mrSges(6,1) - (-t10 * t34 - t15 * t33) * mrSges(6,2) - m(7) * (-t10 * t20 + t15 * t21) - t79 + t91 * (-t10 * t41 + t15 * t44)) * g(1), t93 * (-g(1) * t9 - g(2) * t7 - g(3) * t17) -g(1) * t79 - g(2) * t80 - g(3) * t78];
taug  = t1(:);
