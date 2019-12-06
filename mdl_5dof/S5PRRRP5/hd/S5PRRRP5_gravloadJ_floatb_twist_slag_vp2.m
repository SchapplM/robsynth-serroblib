% Calculate Gravitation load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:47:48
% DurationCPUTime: 0.46s
% Computational Cost: add. (207->64), mult. (311->80), div. (0->0), fcn. (278->8), ass. (0->33)
t69 = mrSges(5,1) + mrSges(6,1);
t67 = mrSges(5,2) + mrSges(6,2);
t22 = qJ(3) + qJ(4);
t19 = cos(t22);
t27 = cos(qJ(3));
t20 = t27 * pkin(3);
t15 = pkin(4) * t19 + t20;
t18 = sin(t22);
t25 = sin(qJ(3));
t66 = mrSges(3,1) + m(6) * (pkin(2) + t15) + m(5) * (t20 + pkin(2)) + m(4) * pkin(2) + t27 * mrSges(4,1) - t25 * mrSges(4,2) + t69 * t19 - t67 * t18;
t29 = -pkin(7) - pkin(6);
t65 = mrSges(3,2) + m(6) * (-qJ(5) + t29) - mrSges(6,3) + m(5) * t29 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t64 = g(1) * t24 + g(2) * t23;
t63 = -m(5) * pkin(3) - mrSges(4,1);
t62 = mrSges(5,1) * t18 + t67 * t19;
t28 = cos(qJ(2));
t46 = t24 * t28;
t11 = -t18 * t46 + t23 * t19;
t43 = t28 * t19;
t61 = -t67 * (-t23 * t18 - t24 * t43) - t69 * t11;
t47 = t23 * t28;
t9 = -t18 * t47 - t24 * t19;
t60 = -t69 * t9 - t67 * (t24 * t18 - t23 * t43);
t57 = m(6) * pkin(4);
t26 = sin(qJ(2));
t50 = g(3) * t26;
t49 = t25 * pkin(3);
t45 = t25 * t28;
t44 = t27 * t28;
t14 = -pkin(4) * t18 - t49;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (-g(3) * t66 + t64 * t65) * t28 + (t65 * g(3) + t64 * t66) * t26, (m(5) * t49 - m(6) * t14 + mrSges(4,1) * t25 + mrSges(6,1) * t18 + mrSges(4,2) * t27 + t62) * t50 + (-(-t23 * t44 + t24 * t25) * mrSges(4,2) - m(6) * (t14 * t47 - t24 * t15) + t63 * (-t23 * t45 - t24 * t27) + t60) * g(2) + (-(-t23 * t25 - t24 * t44) * mrSges(4,2) - m(6) * (t14 * t46 + t23 * t15) + t63 * (t23 * t27 - t24 * t45) + t61) * g(1), (-(-mrSges(6,1) - t57) * t18 + t62) * t50 + (-t9 * t57 + t60) * g(2) + (-t11 * t57 + t61) * g(1), (g(3) * t28 - t64 * t26) * m(6)];
taug = t1(:);
