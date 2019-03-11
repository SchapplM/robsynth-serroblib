% Calculate Gravitation load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:44
% EndTime: 2019-03-08 19:24:45
% DurationCPUTime: 0.84s
% Computational Cost: add. (548->96), mult. (1192->148), div. (0->0), fcn. (1438->14), ass. (0->50)
t44 = sin(qJ(6));
t47 = cos(qJ(6));
t96 = -m(7) * pkin(5) - t47 * mrSges(7,1) + t44 * mrSges(7,2) - mrSges(6,1);
t94 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t97 = m(6) + m(7);
t98 = -m(4) - m(5);
t69 = t97 - t98;
t38 = sin(pkin(11));
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t75 = cos(pkin(11));
t24 = -t49 * t38 - t46 * t75;
t39 = sin(pkin(10));
t41 = cos(pkin(10));
t42 = cos(pkin(6));
t79 = t42 * t49;
t101 = -t39 * t46 + t41 * t79;
t57 = -t46 * t38 + t49 * t75;
t76 = t24 * t42;
t14 = t39 * t76 + t41 * t57;
t45 = sin(qJ(4));
t40 = sin(pkin(6));
t48 = cos(qJ(4));
t84 = t40 * t48;
t100 = -t14 * t45 + t39 * t84;
t22 = t24 * t40;
t99 = t22 * t45 + t42 * t48;
t9 = -t39 * t57 + t41 * t76;
t59 = -t41 * t84 + t45 * t9;
t37 = qJ(4) + pkin(12);
t35 = sin(t37);
t36 = cos(t37);
t95 = -m(5) * pkin(3) - t48 * mrSges(5,1) + t45 * mrSges(5,2) + t94 * t35 + t96 * t36 - mrSges(4,1);
t93 = m(5) * pkin(8) + t44 * mrSges(7,1) + t47 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t88 = t39 * t40;
t86 = t40 * t41;
t85 = t40 * t45;
t81 = t42 * t46;
t31 = t40 * t49 * pkin(2);
t65 = t101 * pkin(2);
t52 = t42 * t57;
t43 = -qJ(5) - pkin(8);
t34 = pkin(4) * t48 + pkin(3);
t21 = t57 * t40;
t16 = -t22 * t36 + t35 * t42;
t13 = t24 * t41 - t39 * t52;
t10 = t39 * t24 + t41 * t52;
t4 = t14 * t36 + t35 * t88;
t2 = -t35 * t86 - t36 * t9;
t1 = [(-m(2) - m(3) - t69) * g(3) (-(mrSges(3,1) * t49 - mrSges(3,2) * t46) * t40 - t97 * (t21 * t34 + t22 * t43 + t31) + t98 * t31 + t93 * t22 + t95 * t21) * g(3) + (-t101 * mrSges(3,1) - (-t39 * t49 - t41 * t81) * mrSges(3,2) + t98 * t65 - t97 * (t10 * t34 + t9 * t43 + t65) + t93 * t9 + t95 * t10) * g(2) + (-(t39 * t81 - t41 * t49) * mrSges(3,2) - t97 * (t13 * t34 - t14 * t43) + t95 * t13 - t93 * t14 + (-t69 * pkin(2) - mrSges(3,1)) * (-t39 * t79 - t41 * t46)) * g(1) (-g(3) * t42 + (-g(1) * t39 + g(2) * t41) * t40) * t69 (-t99 * mrSges(5,1) - (t22 * t48 - t42 * t45) * mrSges(5,2) + t94 * t16 + t96 * (t22 * t35 + t36 * t42)) * g(3) + (-t59 * mrSges(5,1) - (t41 * t85 + t48 * t9) * mrSges(5,2) + t96 * (t35 * t9 - t36 * t86) + t94 * t2) * g(2) + (-t100 * mrSges(5,1) - (-t14 * t48 - t39 * t85) * mrSges(5,2) + t94 * t4 + t96 * (-t14 * t35 + t36 * t88)) * g(1) + (-g(1) * t100 - g(2) * t59 - g(3) * t99) * t97 * pkin(4), t97 * (g(1) * t13 + g(2) * t10 + g(3) * t21) -g(1) * ((-t13 * t47 - t4 * t44) * mrSges(7,1) + (t13 * t44 - t4 * t47) * mrSges(7,2)) - g(2) * ((-t10 * t47 - t2 * t44) * mrSges(7,1) + (t10 * t44 - t2 * t47) * mrSges(7,2)) - g(3) * ((-t16 * t44 - t21 * t47) * mrSges(7,1) + (-t16 * t47 + t21 * t44) * mrSges(7,2))];
taug  = t1(:);
