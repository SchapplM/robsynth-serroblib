% Calculate Gravitation load on the joints for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:12
% EndTime: 2019-03-09 02:00:14
% DurationCPUTime: 0.54s
% Computational Cost: add. (401->77), mult. (360->88), div. (0->0), fcn. (324->10), ass. (0->46)
t67 = mrSges(6,1) + mrSges(7,1);
t66 = -mrSges(6,2) + mrSges(7,3);
t65 = -mrSges(6,3) - mrSges(7,2);
t25 = sin(qJ(5));
t27 = cos(qJ(5));
t64 = t66 * t25 + t27 * t67;
t23 = cos(pkin(10));
t63 = -mrSges(3,1) - m(4) * pkin(2) - t23 * mrSges(4,1) + sin(pkin(10)) * mrSges(4,2);
t62 = m(3) + m(4);
t61 = -m(6) - m(7);
t20 = pkin(10) + qJ(4);
t15 = sin(t20);
t17 = cos(t20);
t60 = t17 * pkin(4) + t15 * pkin(8);
t59 = -m(5) + t61;
t34 = pkin(5) * t27 + qJ(6) * t25;
t58 = (mrSges(5,2) + t65) * t17 + (-m(7) * (-pkin(4) - t34) + m(6) * pkin(4) + mrSges(5,1) + t64) * t15;
t38 = t17 * mrSges(5,1) - t15 * mrSges(5,2);
t57 = t65 * t15 - t38;
t56 = pkin(8) * t61;
t55 = m(7) * pkin(5) + t67;
t54 = m(7) * qJ(6) + t66;
t53 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t52 = g(3) * t15;
t26 = sin(qJ(1));
t51 = t26 * pkin(1);
t28 = cos(qJ(1));
t19 = t28 * pkin(1);
t21 = qJ(1) + pkin(9);
t18 = cos(t21);
t48 = t15 * t18;
t16 = sin(t21);
t47 = t16 * t25;
t46 = t16 * t27;
t45 = t17 * t18;
t44 = t18 * t25;
t43 = t18 * t27;
t41 = m(4) - t59;
t14 = t23 * pkin(3) + pkin(2);
t24 = -pkin(7) - qJ(3);
t40 = t18 * t14 - t16 * t24 + t19;
t4 = t17 * t43 + t47;
t3 = t17 * t44 - t46;
t2 = t17 * t46 - t44;
t1 = t17 * t47 + t43;
t5 = [(-m(5) * t40 - t28 * mrSges(2,1) + t26 * mrSges(2,2) + t65 * t48 + t61 * (pkin(4) * t45 + pkin(8) * t48 + t40) - t55 * t4 - t54 * t3 - t62 * t19 + (-t38 + t63) * t18 + t53 * t16) * g(2) + (t26 * mrSges(2,1) + t28 * mrSges(2,2) + t62 * t51 + t59 * (-t18 * t24 - t51) + t55 * t2 + t54 * t1 + t53 * t18 + (m(5) * t14 + t61 * (-t14 - t60) - t57 - t63) * t16) * g(1) (-m(3) - t41) * g(3) (-g(1) * t16 + g(2) * t18) * t41 (t61 * t60 + (-m(7) * t34 - t64) * t17 + t57) * g(3) + (t58 * t18 + t45 * t56) * g(1) + (t17 * t56 + t58) * g(2) * t16 (t55 * t25 - t54 * t27) * t52 + (t55 * t1 - t54 * t2) * g(2) + (t55 * t3 - t54 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t25 * t52) * m(7)];
taug  = t5(:);
