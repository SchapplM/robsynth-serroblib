% Calculate Gravitation load on the joints for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:26
% EndTime: 2019-12-05 15:53:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (189->51), mult. (246->64), div. (0->0), fcn. (211->10), ass. (0->25)
t14 = pkin(9) + qJ(4);
t10 = cos(t14);
t17 = cos(pkin(9));
t11 = qJ(5) + t14;
t6 = sin(t11);
t7 = cos(t11);
t8 = t17 * pkin(3) + pkin(2);
t9 = sin(t14);
t48 = mrSges(3,1) + m(4) * pkin(2) + t17 * mrSges(4,1) - sin(pkin(9)) * mrSges(4,2) + m(5) * t8 + t10 * mrSges(5,1) - t9 * mrSges(5,2) + m(6) * (pkin(4) * t10 + t8) + t7 * mrSges(6,1) - t6 * mrSges(6,2);
t19 = -pkin(6) - qJ(3);
t47 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(6) * (-pkin(7) + t19) - mrSges(6,3) + m(5) * t19 - mrSges(5,3);
t16 = sin(pkin(8));
t18 = cos(pkin(8));
t46 = g(1) * t18 + g(2) * t16;
t45 = m(6) * pkin(4) + mrSges(5,1);
t21 = cos(qJ(2));
t37 = t16 * t21;
t42 = (-t18 * t7 - t6 * t37) * mrSges(6,1) + (t18 * t6 - t7 * t37) * mrSges(6,2);
t36 = t18 * t21;
t41 = (t16 * t7 - t6 * t36) * mrSges(6,1) + (-t16 * t6 - t7 * t36) * mrSges(6,2);
t20 = sin(qJ(2));
t38 = g(3) * t20;
t35 = m(4) + m(5) + m(6);
t31 = -mrSges(6,1) * t6 - mrSges(6,2) * t7;
t1 = [(-m(2) - m(3) - t35) * g(3), (-t48 * g(3) + t46 * t47) * t21 + (t47 * g(3) + t46 * t48) * t20, (t21 * g(3) - t46 * t20) * t35, (mrSges(5,2) * t10 + t45 * t9 - t31) * t38 + (-(-t10 * t37 + t18 * t9) * mrSges(5,2) - t42 - t45 * (-t18 * t10 - t9 * t37)) * g(2) + (-(-t10 * t36 - t16 * t9) * mrSges(5,2) - t41 - t45 * (t16 * t10 - t9 * t36)) * g(1), -g(1) * t41 - g(2) * t42 - t31 * t38];
taug = t1(:);
