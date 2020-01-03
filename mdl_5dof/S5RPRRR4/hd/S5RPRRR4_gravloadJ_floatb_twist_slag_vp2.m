% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:03
% DurationCPUTime: 0.19s
% Computational Cost: add. (265->52), mult. (143->54), div. (0->0), fcn. (102->10), ass. (0->30)
t48 = mrSges(5,2) - mrSges(6,3);
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t47 = mrSges(6,1) * t30 - mrSges(6,2) * t28;
t46 = -mrSges(5,1) - t47;
t27 = qJ(1) + pkin(9);
t24 = qJ(3) + t27;
t21 = qJ(4) + t24;
t15 = sin(t21);
t16 = cos(t21);
t45 = t16 * pkin(4) + t15 * pkin(8);
t19 = sin(t24);
t13 = pkin(3) * t19;
t20 = cos(t24);
t14 = pkin(3) * t20;
t22 = sin(t27);
t29 = sin(qJ(1));
t42 = t29 * pkin(1) + pkin(2) * t22;
t23 = cos(t27);
t31 = cos(qJ(1));
t41 = t31 * pkin(1) + pkin(2) * t23;
t40 = t14 + t41;
t39 = -m(3) * pkin(1) - mrSges(2,1);
t38 = t15 * pkin(4) - pkin(8) * t16;
t37 = t13 + t38;
t35 = t48 * t15 + t46 * t16;
t34 = t46 * t15 - t48 * t16;
t33 = -t19 * mrSges(4,1) - t20 * mrSges(4,2) + t34;
t32 = -t20 * mrSges(4,1) + t19 * mrSges(4,2) + t35;
t1 = [(-mrSges(2,2) * t31 - mrSges(3,1) * t22 - mrSges(3,2) * t23 - m(4) * t42 - m(5) * (t13 + t42) - m(6) * (t37 + t42) + t39 * t29 + t33) * g(3) + (t29 * mrSges(2,2) - t23 * mrSges(3,1) + t22 * mrSges(3,2) - m(4) * t41 - m(5) * t40 - m(6) * (t40 + t45) + t39 * t31 + t32) * g(2), (-m(3) - m(4) - m(5) - m(6)) * g(1), (-m(5) * t13 - m(6) * t37 + t33) * g(3) + (-m(5) * t14 - m(6) * (t14 + t45) + t32) * g(2), (-m(6) * t38 + t34) * g(3) + (-m(6) * t45 + t35) * g(2), -g(1) * t47 + (g(2) * t15 - g(3) * t16) * (mrSges(6,1) * t28 + mrSges(6,2) * t30)];
taug = t1(:);
