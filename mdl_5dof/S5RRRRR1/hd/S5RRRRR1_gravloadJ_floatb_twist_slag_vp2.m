% Calculate Gravitation load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:55
% EndTime: 2019-03-08 18:36:57
% DurationCPUTime: 0.42s
% Computational Cost: add. (359->81), mult. (350->95), div. (0->0), fcn. (292->10), ass. (0->52)
t23 = qJ(2) + qJ(3);
t21 = qJ(4) + t23;
t16 = sin(t21);
t17 = cos(t21);
t24 = sin(qJ(5));
t65 = mrSges(6,2) * t24;
t83 = m(6) * pkin(6) + mrSges(6,3);
t84 = t16 * t65 + t17 * t83;
t19 = sin(t23);
t20 = cos(t23);
t66 = mrSges(5,2) * t17;
t82 = mrSges(4,1) * t19 + mrSges(5,1) * t16 + mrSges(4,2) * t20 + t66;
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t81 = g(1) * t29 + g(2) * t26;
t14 = t16 * mrSges(5,2);
t80 = t20 * mrSges(4,1) + t17 * mrSges(5,1) - t19 * mrSges(4,2) + t16 * mrSges(6,3) - t14;
t10 = t17 * t65;
t27 = cos(qJ(5));
t57 = t27 * mrSges(6,1);
t76 = t17 * t57 - t10 + t80;
t74 = mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t28 = cos(qJ(2));
t22 = t28 * pkin(2);
t25 = sin(qJ(2));
t46 = -t28 * mrSges(3,1) + t25 * mrSges(3,2);
t47 = t17 * pkin(4) + t16 * pkin(6);
t71 = pkin(3) * t20;
t52 = t22 + t71;
t7 = pkin(1) + t52;
t73 = m(6) * (t47 + t7) + mrSges(2,1) + m(4) * (t22 + pkin(1)) + m(5) * t7 + m(3) * pkin(1) - t46 + t80;
t72 = pkin(3) * t19;
t70 = pkin(4) * t16;
t67 = t25 * pkin(2);
t59 = t26 * t24;
t58 = t26 * t27;
t56 = t29 * t24;
t55 = t29 * t27;
t53 = t16 * t57;
t49 = t84 * t26;
t48 = t84 * t29;
t41 = m(6) * pkin(4) + mrSges(5,1) + t57;
t39 = -t47 - t71;
t8 = -t67 - t72;
t32 = m(6) * (t8 - t70) - t53;
t31 = m(6) * (-t70 - t72) - t53;
t30 = t41 * t16 + t66;
t4 = t17 * t55 - t59;
t3 = -t17 * t56 - t58;
t2 = -t17 * t58 - t56;
t1 = t17 * t59 - t55;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) + t74 * t26 - t73 * t29) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + t73 * t26 + t74 * t29) * g(1), -g(1) * (t32 * t29 + t48) - g(2) * (t32 * t26 + t49) + (-t46 + m(4) * t22 + m(5) * t52 - m(6) * (-t22 + t39) + t76) * g(3) + t81 * (m(4) * t67 - m(5) * t8 + mrSges(3,1) * t25 + mrSges(3,2) * t28 + t82) -g(1) * (t31 * t29 + t48) - g(2) * (t31 * t26 + t49) + (m(5) * t71 - m(6) * t39 + t76) * g(3) + (m(5) * t72 + t82) * t81 (t83 * t16 + t41 * t17 - t10 - t14) * g(3) + (t30 * t26 - t49) * g(2) + (t30 * t29 - t48) * g(1), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (mrSges(6,1) * t24 + mrSges(6,2) * t27) * t16];
taug  = t5(:);
