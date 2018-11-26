% Calculate Gravitation load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2018-11-23 14:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:57:52
% EndTime: 2018-11-23 14:57:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (926->88), mult. (942->119), div. (0->0), fcn. (881->16), ass. (0->52)
t90 = m(6) + m(7);
t95 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t47 = sin(qJ(6));
t49 = cos(qJ(6));
t94 = mrSges(7,1) * t47 + mrSges(7,2) * t49 - mrSges(5,2) + mrSges(6,3);
t91 = m(7) * pkin(9) + pkin(4) * t90 + t95;
t87 = -t90 * qJ(5) - t94;
t42 = pkin(11) + qJ(4);
t40 = sin(t42);
t41 = cos(t42);
t45 = cos(pkin(11));
t86 = m(4) * pkin(2) + t45 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + mrSges(3,1) + t95 * t41 + t94 * t40;
t85 = m(4) * qJ(3) + m(7) * pkin(5) + t49 * mrSges(7,1) - t47 * mrSges(7,2) + mrSges(6,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t44 = sin(pkin(10));
t48 = sin(qJ(2));
t70 = pkin(6) + qJ(2);
t61 = cos(t70) / 0.2e1;
t71 = pkin(6) - qJ(2);
t63 = cos(t71);
t51 = t63 / 0.2e1 + t61;
t73 = cos(pkin(10));
t18 = t44 * t48 - t73 * t51;
t39 = pkin(3) * t45 + pkin(2);
t46 = -pkin(8) - qJ(3);
t60 = sin(t70) / 0.2e1;
t62 = sin(t71);
t31 = t60 - t62 / 0.2e1;
t50 = cos(qJ(2));
t56 = t73 * t31 + t44 * t50;
t84 = -t18 * t39 - t46 * t56;
t83 = t18 * t41;
t21 = t44 * t51 + t73 * t48;
t82 = t21 * t41;
t30 = t60 + t62 / 0.2e1;
t81 = t30 * t41;
t55 = -t44 * t31 + t73 * t50;
t77 = -t21 * t39 - t46 * t55;
t32 = t61 - t63 / 0.2e1;
t76 = t30 * t39 + t32 * t46;
t75 = qJ(5) * t40;
t74 = cos(pkin(6));
t72 = sin(pkin(6));
t69 = m(4) + m(5) + t90;
t68 = -pkin(4) * t83 - t18 * t75 + t84;
t66 = t44 * t72;
t65 = -pkin(4) * t82 - t21 * t75 + t77;
t64 = pkin(4) * t81 + t30 * t75 + t76;
t57 = t73 * t72;
t14 = -t32 * t40 - t74 * t41;
t5 = t40 * t55 - t41 * t66;
t3 = t40 * t56 + t41 * t57;
t1 = [(-m(2) - m(3) - t69) * g(3) (-m(5) * t76 - m(6) * t64 - m(7) * (pkin(9) * t81 + t64) + t85 * t32 - t86 * t30) * g(3) + (-m(5) * t84 - m(6) * t68 - m(7) * (-pkin(9) * t83 + t68) - t85 * t56 + t86 * t18) * g(2) + (-m(5) * t77 - m(6) * t65 - m(7) * (-pkin(9) * t82 + t65) - t85 * t55 + t86 * t21) * g(1) (-g(1) * t21 - g(2) * t18 + g(3) * t30) * t69 (t87 * (-t32 * t41 + t74 * t40) + t91 * t14) * g(3) + (t87 * (-t40 * t57 + t41 * t56) + t91 * t3) * g(2) + (t87 * (t40 * t66 + t41 * t55) + t91 * t5) * g(1), t90 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((-t21 * t47 + t49 * t5) * mrSges(7,1) + (-t21 * t49 - t47 * t5) * mrSges(7,2)) - g(2) * ((-t18 * t47 + t3 * t49) * mrSges(7,1) + (-t18 * t49 - t3 * t47) * mrSges(7,2)) - g(3) * ((t14 * t49 + t30 * t47) * mrSges(7,1) + (-t14 * t47 + t30 * t49) * mrSges(7,2))];
taug  = t1(:);
