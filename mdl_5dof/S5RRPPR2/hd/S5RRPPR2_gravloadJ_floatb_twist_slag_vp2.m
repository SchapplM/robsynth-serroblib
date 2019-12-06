% Calculate Gravitation load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:03
% EndTime: 2019-12-05 18:20:04
% DurationCPUTime: 0.30s
% Computational Cost: add. (280->61), mult. (198->71), div. (0->0), fcn. (168->10), ass. (0->34)
t55 = mrSges(4,2) - mrSges(5,3);
t21 = sin(pkin(9));
t22 = cos(pkin(9));
t54 = mrSges(4,1) - m(6) * (-pkin(4) * t22 - pkin(3)) + t22 * mrSges(5,1) + (m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3)) * t21;
t52 = m(5) + m(6);
t20 = qJ(1) + qJ(2);
t17 = pkin(8) + t20;
t15 = sin(t17);
t16 = cos(t17);
t18 = sin(t20);
t19 = cos(t20);
t25 = cos(qJ(5));
t23 = sin(qJ(5));
t43 = t22 * t23;
t7 = -t15 * t25 + t16 * t43;
t42 = t22 * t25;
t8 = -t15 * t23 - t16 * t42;
t51 = t19 * mrSges(3,1) - t8 * mrSges(6,1) - t18 * mrSges(3,2) - t7 * mrSges(6,2) - t55 * t15 + t54 * t16;
t5 = t15 * t43 + t16 * t25;
t6 = t15 * t42 - t16 * t23;
t50 = t18 * mrSges(3,1) + t6 * mrSges(6,1) + t19 * mrSges(3,2) - t5 * mrSges(6,2) + t54 * t15 + t55 * t16;
t49 = pkin(2) * t18;
t48 = pkin(2) * t19;
t24 = sin(qJ(1));
t47 = t24 * pkin(1);
t26 = cos(qJ(1));
t46 = t26 * pkin(1);
t11 = t16 * qJ(4);
t41 = t11 - t49;
t38 = -t47 - t49;
t34 = -t15 * qJ(4) - t48;
t32 = -t15 * pkin(3) + t41;
t30 = -t16 * pkin(3) + t34;
t1 = [(t24 * mrSges(2,1) + t26 * mrSges(2,2) + m(3) * t47 - m(4) * t38 - m(5) * (t32 - t47) - m(6) * (t11 + t38) + t50) * g(3) + (t26 * mrSges(2,1) - t24 * mrSges(2,2) + m(3) * t46 - m(4) * (-t46 - t48) - m(5) * (t30 - t46) - m(6) * (t34 - t46) + t51) * g(2), (m(4) * t49 - m(5) * t32 - m(6) * t41 + t50) * g(3) + (m(4) * t48 - m(5) * t30 - m(6) * t34 + t51) * g(2), (-m(4) - t52) * g(1), t52 * (-g(2) * t16 - g(3) * t15), -g(2) * (t5 * mrSges(6,1) + t6 * mrSges(6,2)) - g(3) * (-t7 * mrSges(6,1) + t8 * mrSges(6,2)) - g(1) * (-mrSges(6,1) * t23 - mrSges(6,2) * t25) * t21];
taug = t1(:);
