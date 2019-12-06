% Calculate Gravitation load on the joints for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:26
% EndTime: 2019-12-05 16:02:28
% DurationCPUTime: 0.51s
% Computational Cost: add. (207->56), mult. (518->87), div. (0->0), fcn. (582->10), ass. (0->33)
t52 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t60 = m(5) + m(6);
t51 = m(4) + t60;
t23 = sin(qJ(5));
t26 = cos(qJ(5));
t54 = -m(6) * pkin(4) - t26 * mrSges(6,1) + t23 * mrSges(6,2) - mrSges(5,1);
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t61 = t54 * t24 - t52 * t27 + mrSges(3,2) - mrSges(4,3);
t58 = -t23 * mrSges(6,1) - t26 * mrSges(6,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t57 = pkin(2) * t51 + t60 * pkin(7) - t58;
t47 = -t51 * qJ(3) + t61;
t21 = sin(pkin(5));
t46 = t21 * t24;
t25 = sin(qJ(2));
t45 = t21 * t25;
t44 = t21 * t27;
t28 = cos(qJ(2));
t43 = t21 * t28;
t42 = pkin(2) * t43 + qJ(3) * t45;
t41 = cos(pkin(5));
t38 = t25 * t41;
t37 = t28 * t41;
t22 = cos(pkin(9));
t20 = sin(pkin(9));
t13 = -t24 * t43 + t41 * t27;
t11 = -t20 * t38 + t22 * t28;
t10 = t20 * t37 + t22 * t25;
t9 = t20 * t28 + t22 * t38;
t8 = t20 * t25 - t22 * t37;
t4 = t22 * t44 - t8 * t24;
t2 = t10 * t24 + t20 * t44;
t1 = [(-m(2) - m(3) - t51) * g(3), (-m(4) * t42 - t60 * (pkin(7) * t43 + t42) + (t61 * t25 + t58 * t28) * t21) * g(3) + (t47 * t9 + t57 * t8) * g(2) + (t57 * t10 + t47 * t11) * g(1), t51 * (-g(1) * t10 - g(2) * t8 + g(3) * t43), (t52 * t13 + t54 * (-t41 * t24 - t27 * t43)) * g(3) + (-t52 * t4 + t54 * (t22 * t46 + t8 * t27)) * g(2) + (t52 * t2 + t54 * (t10 * t27 - t20 * t46)) * g(1), -g(1) * ((t11 * t26 - t2 * t23) * mrSges(6,1) + (-t11 * t23 - t2 * t26) * mrSges(6,2)) - g(2) * ((t4 * t23 + t9 * t26) * mrSges(6,1) + (-t9 * t23 + t4 * t26) * mrSges(6,2)) - g(3) * ((-t13 * t23 + t26 * t45) * mrSges(6,1) + (-t13 * t26 - t23 * t45) * mrSges(6,2))];
taug = t1(:);
