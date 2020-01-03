% Calculate Gravitation load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:22
% DurationCPUTime: 0.20s
% Computational Cost: add. (284->55), mult. (156->56), div. (0->0), fcn. (112->10), ass. (0->36)
t53 = mrSges(5,2) - mrSges(6,3);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t52 = mrSges(6,1) * t33 - mrSges(6,2) * t31;
t51 = -mrSges(5,1) - t52;
t30 = qJ(1) + qJ(2);
t25 = pkin(9) + t30;
t24 = qJ(4) + t25;
t15 = sin(t24);
t16 = cos(t24);
t50 = t16 * pkin(4) + t15 * pkin(8);
t26 = sin(t30);
t22 = pkin(2) * t26;
t27 = cos(t30);
t23 = pkin(2) * t27;
t20 = sin(t25);
t13 = pkin(3) * t20;
t47 = t13 + t22;
t21 = cos(t25);
t14 = pkin(3) * t21;
t46 = t14 + t23;
t32 = sin(qJ(1));
t28 = t32 * pkin(1);
t45 = t22 + t28;
t34 = cos(qJ(1));
t29 = t34 * pkin(1);
t44 = t23 + t29;
t43 = -m(3) * pkin(1) - mrSges(2,1);
t42 = t15 * pkin(4) - pkin(8) * t16;
t41 = t46 + t50;
t39 = t42 + t47;
t38 = t53 * t15 + t51 * t16;
t37 = t51 * t15 - t53 * t16;
t36 = -t26 * mrSges(3,1) - t20 * mrSges(4,1) - t27 * mrSges(3,2) - t21 * mrSges(4,2) + t37;
t35 = -t27 * mrSges(3,1) - t21 * mrSges(4,1) + t26 * mrSges(3,2) + t20 * mrSges(4,2) + t38;
t1 = [(-mrSges(2,2) * t34 - m(4) * t45 - m(5) * (t13 + t45) - m(6) * (t28 + t39) + t43 * t32 + t36) * g(3) + (t32 * mrSges(2,2) - m(4) * t44 - m(5) * (t14 + t44) - m(6) * (t29 + t41) + t43 * t34 + t35) * g(2), (-m(4) * t22 - m(5) * t47 - m(6) * t39 + t36) * g(3) + (-m(4) * t23 - m(5) * t46 - m(6) * t41 + t35) * g(2), (-m(4) - m(5) - m(6)) * g(1), (-m(6) * t42 + t37) * g(3) + (-m(6) * t50 + t38) * g(2), -g(1) * t52 + (g(2) * t15 - g(3) * t16) * (mrSges(6,1) * t31 + mrSges(6,2) * t33)];
taug = t1(:);
