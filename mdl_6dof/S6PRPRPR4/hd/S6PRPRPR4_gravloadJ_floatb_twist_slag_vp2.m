% Calculate Gravitation load on the joints for
% S6PRPRPR4
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:57:06
% EndTime: 2018-11-23 14:57:07
% DurationCPUTime: 0.67s
% Computational Cost: add. (982->81), mult. (976->106), div. (0->0), fcn. (920->18), ass. (0->47)
t36 = pkin(12) + qJ(6);
t32 = sin(t36);
t34 = cos(t36);
t38 = sin(pkin(12));
t41 = cos(pkin(12));
t80 = mrSges(5,1) + m(7) * (pkin(5) * t41 + pkin(4)) + t34 * mrSges(7,1) - t32 * mrSges(7,2) + m(6) * pkin(4) + t41 * mrSges(6,1) - t38 * mrSges(6,2);
t79 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t37 = pkin(11) + qJ(4);
t33 = sin(t37);
t35 = cos(t37);
t42 = cos(pkin(11));
t78 = m(4) * pkin(2) + t42 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + mrSges(3,1) + t80 * t35 - t79 * t33;
t85 = m(6) + m(7);
t82 = -m(5) - t85;
t77 = -m(4) * qJ(3) - t32 * mrSges(7,1) - t41 * mrSges(6,2) - t34 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t38;
t73 = cos(pkin(6));
t72 = cos(pkin(10));
t71 = sin(pkin(6));
t70 = pkin(6) - qJ(2);
t69 = pkin(6) + qJ(2);
t68 = m(4) - t82;
t40 = sin(pkin(10));
t64 = t40 * t71;
t63 = cos(t70);
t62 = sin(t70);
t61 = cos(t69) / 0.2e1;
t60 = sin(t69) / 0.2e1;
t57 = t72 * t71;
t22 = t60 - t62 / 0.2e1;
t46 = cos(qJ(2));
t55 = t72 * t22 + t40 * t46;
t54 = -t40 * t22 + t72 * t46;
t49 = t63 / 0.2e1 + t61;
t45 = sin(qJ(2));
t44 = -pkin(8) - qJ(3);
t31 = pkin(3) * t42 + pkin(2);
t23 = t61 - t63 / 0.2e1;
t21 = t60 + t62 / 0.2e1;
t14 = t40 * t49 + t72 * t45;
t11 = t40 * t45 - t72 * t49;
t8 = -t23 * t35 + t73 * t33;
t7 = -t23 * t33 - t73 * t35;
t4 = t33 * t64 + t35 * t54;
t3 = t33 * t54 - t35 * t64;
t2 = -t33 * t57 + t35 * t55;
t1 = t33 * t55 + t35 * t57;
t5 = [(-m(2) - m(3) - t68) * g(3) (t82 * (t21 * t31 + t23 * t44) - t77 * t23 - t78 * t21) * g(3) + (t82 * (-t11 * t31 - t44 * t55) + t77 * t55 + t78 * t11) * g(2) + (t82 * (-t14 * t31 - t44 * t54) + t77 * t54 + t78 * t14) * g(1) (-g(1) * t14 - g(2) * t11 + g(3) * t21) * t68 (t80 * t7 + t79 * t8) * g(3) + (t80 * t1 + t79 * t2) * g(2) + (t80 * t3 + t79 * t4) * g(1), t85 * (-g(1) * t3 - g(2) * t1 - g(3) * t7) -g(1) * ((t14 * t34 - t32 * t4) * mrSges(7,1) + (-t14 * t32 - t34 * t4) * mrSges(7,2)) - g(2) * ((t11 * t34 - t2 * t32) * mrSges(7,1) + (-t11 * t32 - t2 * t34) * mrSges(7,2)) - g(3) * ((-t21 * t34 - t32 * t8) * mrSges(7,1) + (t21 * t32 - t34 * t8) * mrSges(7,2))];
taug  = t5(:);
