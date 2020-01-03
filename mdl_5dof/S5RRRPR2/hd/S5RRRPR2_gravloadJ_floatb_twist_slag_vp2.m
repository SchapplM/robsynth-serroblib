% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:22
% DurationCPUTime: 0.21s
% Computational Cost: add. (307->58), mult. (167->58), div. (0->0), fcn. (120->10), ass. (0->34)
t51 = mrSges(5,2) - mrSges(6,3);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t50 = mrSges(6,1) * t33 - mrSges(6,2) * t31;
t49 = -mrSges(5,1) - t50;
t30 = qJ(1) + qJ(2);
t25 = sin(t30);
t20 = pkin(2) * t25;
t26 = cos(t30);
t21 = pkin(2) * t26;
t27 = qJ(3) + t30;
t23 = sin(t27);
t15 = pkin(3) * t23;
t24 = cos(t27);
t16 = pkin(3) * t24;
t32 = sin(qJ(1));
t28 = t32 * pkin(1);
t46 = t20 + t28;
t34 = cos(qJ(1));
t29 = t34 * pkin(1);
t45 = t21 + t29;
t22 = pkin(9) + t27;
t13 = sin(t22);
t14 = cos(t22);
t44 = t14 * pkin(4) + t13 * pkin(8) + t16;
t43 = -m(3) * pkin(1) - mrSges(2,1);
t42 = t21 + t44;
t41 = t13 * pkin(4) - pkin(8) * t14 + t15;
t39 = t20 + t41;
t38 = -t23 * mrSges(4,1) - t24 * mrSges(4,2) + t49 * t13 - t51 * t14;
t37 = -t24 * mrSges(4,1) + t23 * mrSges(4,2) + t51 * t13 + t49 * t14;
t36 = -t25 * mrSges(3,1) - t26 * mrSges(3,2) + t38;
t35 = -t26 * mrSges(3,1) + t25 * mrSges(3,2) + t37;
t1 = [(-mrSges(2,2) * t34 - m(4) * t46 - m(5) * (t15 + t46) - m(6) * (t28 + t39) + t43 * t32 + t36) * g(3) + (t32 * mrSges(2,2) - m(4) * t45 - m(5) * (t16 + t45) - m(6) * (t29 + t42) + t43 * t34 + t35) * g(2), (-m(4) * t20 - m(5) * (t15 + t20) - m(6) * t39 + t36) * g(3) + (-m(4) * t21 - m(5) * (t16 + t21) - m(6) * t42 + t35) * g(2), (-m(5) * t15 - m(6) * t41 + t38) * g(3) + (-m(5) * t16 - m(6) * t44 + t37) * g(2), (-m(5) - m(6)) * g(1), -g(1) * t50 + (g(2) * t13 - g(3) * t14) * (mrSges(6,1) * t31 + mrSges(6,2) * t33)];
taug = t1(:);
