% Calculate Gravitation load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:38:17
% DurationCPUTime: 0.27s
% Computational Cost: add. (199->49), mult. (165->49), div. (0->0), fcn. (121->10), ass. (0->30)
t20 = qJ(1) + pkin(8);
t13 = cos(t20);
t51 = t13 * g(3);
t47 = m(5) + m(6);
t46 = -m(3) - m(4) - t47;
t50 = pkin(1) * t46 - mrSges(2,1);
t19 = qJ(3) + pkin(9);
t14 = qJ(5) + t19;
t7 = sin(t14);
t8 = cos(t14);
t49 = mrSges(6,1) * t7 + mrSges(6,2) * t8;
t10 = sin(t19);
t12 = cos(t19);
t22 = sin(qJ(3));
t33 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t48 = -t12 * mrSges(5,1) + mrSges(4,2) * t22 + t10 * mrSges(5,2) - t33;
t24 = cos(qJ(3));
t34 = m(5) * pkin(3) + mrSges(4,1);
t45 = -t22 * t34 + m(6) * (-pkin(3) * t22 - pkin(4) * t10) - mrSges(5,1) * t10 - mrSges(4,2) * t24 - mrSges(5,2) * t12;
t16 = t24 * pkin(3);
t37 = pkin(4) * t12 + t16;
t44 = -mrSges(3,1) - m(4) * pkin(2) - mrSges(4,1) * t24 - m(5) * (t16 + pkin(2)) - m(6) * (pkin(2) + t37) + t48;
t21 = -qJ(4) - pkin(6);
t42 = m(4) * pkin(6) - m(5) * t21 - m(6) * (-pkin(7) + t21) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t11 = sin(t20);
t39 = g(2) * t11;
t35 = t49 * t51;
t25 = cos(qJ(1));
t23 = sin(qJ(1));
t1 = [(-mrSges(2,2) * t25 + t44 * t11 + t42 * t13 + t50 * t23) * g(3) + (t23 * mrSges(2,2) - t42 * t11 + t44 * t13 + t50 * t25) * g(2), t46 * g(1), -t35 + (-m(6) * t37 - t24 * t34 + t48) * g(1) + t45 * t51 + (t49 - t45) * t39, t47 * (g(2) * t13 + g(3) * t11), -g(1) * t33 + t39 * t49 - t35];
taug = t1(:);
