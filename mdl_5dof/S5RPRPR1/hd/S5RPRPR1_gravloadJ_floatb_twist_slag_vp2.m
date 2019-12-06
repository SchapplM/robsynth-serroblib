% Calculate Gravitation load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:04
% DurationCPUTime: 0.27s
% Computational Cost: add. (147->52), mult. (183->56), div. (0->0), fcn. (135->8), ass. (0->28)
t18 = cos(qJ(3));
t14 = qJ(3) + pkin(8);
t9 = qJ(5) + t14;
t5 = sin(t9);
t6 = cos(t9);
t25 = -t5 * mrSges(6,1) - t6 * mrSges(6,2);
t16 = sin(qJ(3));
t33 = pkin(3) * t16;
t7 = sin(t14);
t8 = cos(t14);
t41 = m(6) * (pkin(4) * t7 + t33) - t25 + t7 * mrSges(5,1) + t8 * mrSges(5,2) + mrSges(4,2) * t18;
t40 = m(5) + m(6);
t39 = m(3) + m(4) + t40;
t27 = m(5) * pkin(3) + mrSges(4,1);
t38 = -t27 * t18 - m(6) * (pkin(3) * t18 + pkin(4) * t8) - mrSges(5,1) * t8 + mrSges(4,2) * t16 + mrSges(5,2) * t7;
t37 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t36 = -m(5) * t33 - mrSges(4,1) * t16 + mrSges(2,2) - mrSges(3,3) - t41;
t35 = mrSges(6,1) * t6;
t34 = mrSges(6,2) * t5;
t17 = sin(qJ(1));
t32 = g(1) * t17;
t19 = cos(qJ(1));
t31 = g(2) * t19;
t15 = -qJ(4) - pkin(6);
t13 = -pkin(7) + t15;
t4 = t19 * t34;
t3 = t17 * t35;
t1 = [(-t39 * (t19 * pkin(1) + t17 * qJ(2)) + (-m(4) * pkin(6) + m(5) * t15 + m(6) * t13 - t37) * t19 + t36 * t17) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(6)) - m(5) * (-pkin(1) + t15) - m(6) * (-pkin(1) + t13) + t37) * t17 + (-t39 * qJ(2) + t36) * t19) * g(1), (t31 - t32) * t39, -g(1) * t3 - g(2) * t4 + (t27 * t16 + t41) * g(3) + (t35 - t38) * t31 + (t34 + t38) * t32, t40 * (-g(1) * t19 - g(2) * t17), -g(1) * (-t17 * t34 + t3) - g(2) * (-t19 * t35 + t4) - g(3) * t25];
taug = t1(:);
