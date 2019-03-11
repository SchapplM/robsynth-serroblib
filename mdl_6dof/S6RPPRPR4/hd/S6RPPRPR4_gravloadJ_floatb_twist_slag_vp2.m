% Calculate Gravitation load on the joints for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:16
% EndTime: 2019-03-09 01:46:17
% DurationCPUTime: 0.46s
% Computational Cost: add. (258->76), mult. (409->87), div. (0->0), fcn. (435->10), ass. (0->40)
t59 = m(6) + m(7);
t62 = m(7) * pkin(8) + mrSges(7,3);
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t34 = -t24 * mrSges(5,1) + t22 * mrSges(5,2);
t61 = -m(5) * pkin(3) - mrSges(4,1) + t34;
t60 = -m(4) - m(5);
t58 = mrSges(2,1) + mrSges(3,1);
t57 = mrSges(2,2) - mrSges(3,3);
t55 = t59 - t60;
t19 = qJ(4) + pkin(10);
t14 = cos(t19);
t21 = sin(qJ(6));
t23 = cos(qJ(6));
t29 = m(7) * pkin(5) + t23 * mrSges(7,1) - t21 * mrSges(7,2);
t13 = sin(t19);
t32 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t54 = -t29 * t14 - t32;
t52 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t48 = t13 * pkin(8);
t47 = t24 * pkin(4);
t46 = cos(qJ(1));
t45 = sin(qJ(1));
t44 = t13 * mrSges(7,3);
t43 = t14 * t21;
t42 = t14 * t23;
t41 = t46 * pkin(1) + t45 * qJ(2);
t40 = cos(pkin(9));
t39 = sin(pkin(9));
t38 = t46 * pkin(2) + t41;
t35 = -pkin(1) * t45 + t46 * qJ(2);
t31 = mrSges(7,1) * t21 + mrSges(7,2) * t23;
t27 = -pkin(2) * t45 + t35;
t20 = -qJ(5) - pkin(7);
t12 = pkin(3) + t47;
t8 = t39 * t46 - t40 * t45;
t7 = -t39 * t45 - t40 * t46;
t2 = t8 * t21 - t42 * t7;
t1 = t8 * t23 + t43 * t7;
t3 = [(-m(3) * t41 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t58 * t46 + t57 * t45 + t60 * t38 - t59 * (-t7 * t12 - t8 * t20 + t38) + t52 * t8 + (t32 - m(7) * (-t14 * pkin(5) - t48) + t44 - t61) * t7) * g(2) + (-m(3) * t35 + t57 * t46 + t58 * t45 + t60 * t27 - t59 * (t8 * t12 - t7 * t20 + t27) + (-t62 * t13 + t54 + t61) * t8 + (-t31 + t52) * t7) * g(1) (-t45 * g(1) + t46 * g(2)) * (m(3) + t55) t55 * g(3) (-t34 + m(6) * t47 - m(7) * (-t47 - t48) + t44 - t54) * g(3) + (t7 * g(1) + t8 * g(2)) * (-mrSges(5,2) * t24 + (-mrSges(6,2) + t62) * t14 + (-mrSges(6,1) - t29) * t13 + (-t59 * pkin(4) - mrSges(5,1)) * t22) t59 * (-g(1) * t8 + g(2) * t7) -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * ((-t7 * t23 + t43 * t8) * mrSges(7,1) + (t7 * t21 + t42 * t8) * mrSges(7,2)) - g(3) * t31 * t13];
taug  = t3(:);
