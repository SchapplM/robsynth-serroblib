% Calculate Gravitation load on the joints for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2018-11-23 15:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:57:16
% EndTime: 2018-11-23 15:57:17
% DurationCPUTime: 0.62s
% Computational Cost: add. (422->84), mult. (387->98), div. (0->0), fcn. (347->10), ass. (0->50)
t75 = mrSges(6,1) + mrSges(7,1);
t74 = -mrSges(6,2) + mrSges(7,3);
t73 = m(6) + m(7);
t72 = -mrSges(6,3) - mrSges(7,2);
t24 = sin(qJ(5));
t27 = cos(qJ(5));
t71 = t74 * t24 + t75 * t27;
t22 = qJ(1) + pkin(9);
t16 = sin(t22);
t70 = g(2) * t16;
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t41 = t28 * mrSges(4,1) - t25 * mrSges(4,2);
t69 = -m(4) * pkin(2) - mrSges(3,1) - t41;
t68 = m(3) + m(4);
t21 = qJ(3) + pkin(10);
t15 = sin(t21);
t17 = cos(t21);
t66 = -t17 * pkin(4) - t15 * pkin(8);
t36 = pkin(5) * t27 + qJ(6) * t24;
t58 = pkin(3) * t25;
t65 = t73 * t58 + t72 * t17 + (-m(7) * (-pkin(4) - t36) + m(6) * pkin(4) + t71) * t15;
t45 = m(5) + t73;
t39 = t17 * mrSges(5,1) - t15 * mrSges(5,2);
t63 = -t72 * t15 + t39;
t62 = pkin(8) * t73;
t61 = m(7) * pkin(5) + t75;
t60 = m(7) * qJ(6) + t74;
t59 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t55 = g(3) * t15;
t26 = sin(qJ(1));
t54 = t26 * pkin(1);
t19 = t28 * pkin(3);
t29 = cos(qJ(1));
t20 = t29 * pkin(1);
t18 = cos(t22);
t51 = t15 * t18;
t50 = t16 * t24;
t49 = t16 * t27;
t48 = t17 * t18;
t47 = t18 * t24;
t46 = t18 * t27;
t14 = t19 + pkin(2);
t23 = -qJ(4) - pkin(7);
t43 = t18 * t14 - t16 * t23 + t20;
t4 = t17 * t46 + t50;
t3 = t17 * t47 - t49;
t2 = t17 * t49 - t47;
t1 = t17 * t50 + t46;
t5 = [(-m(5) * t43 - t29 * mrSges(2,1) + t26 * mrSges(2,2) + t72 * t51 - t73 * (pkin(4) * t48 + pkin(8) * t51 + t43) - t61 * t4 - t60 * t3 - t68 * t20 + (-t39 + t69) * t18 + t59 * t16) * g(2) + (t26 * mrSges(2,1) + t29 * mrSges(2,2) + t68 * t54 - t45 * (-t18 * t23 - t54) + t61 * t2 + t60 * t1 + t59 * t18 + (m(5) * t14 - t73 * (-t14 + t66) + t63 - t69) * t16) * g(1) (-t45 - t68) * g(3) (t65 * t18 - t48 * t62) * g(1) + (-m(5) * t19 - t41 - t73 * (t19 - t66) + (-m(7) * t36 - t71) * t17 - t63) * g(3) + (m(5) * t58 + mrSges(4,1) * t25 + mrSges(5,1) * t15 + mrSges(4,2) * t28 + mrSges(5,2) * t17) * (g(1) * t18 + t70) + (-t17 * t62 + t65) * t70 (-g(1) * t16 + g(2) * t18) * t45 (t61 * t24 - t60 * t27) * t55 + (t61 * t1 - t60 * t2) * g(2) + (t61 * t3 - t60 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t24 * t55) * m(7)];
taug  = t5(:);
