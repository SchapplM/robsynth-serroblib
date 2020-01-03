% Calculate Gravitation load on the joints for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:09
% EndTime: 2019-12-31 18:42:11
% DurationCPUTime: 0.44s
% Computational Cost: add. (235->64), mult. (281->78), div. (0->0), fcn. (248->8), ass. (0->36)
t49 = -mrSges(5,1) - mrSges(6,1);
t47 = mrSges(5,2) + mrSges(6,2);
t53 = mrSges(5,3) + mrSges(6,3);
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t8 = t17 * pkin(4) + pkin(3);
t52 = -m(5) * pkin(3) - m(6) * t8 + t47 * t14 + t49 * t17;
t13 = -qJ(5) - pkin(7);
t51 = -m(5) * pkin(7) + m(6) * t13 - t53;
t12 = qJ(1) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t50 = g(1) * t10 + g(2) * t9;
t42 = m(6) * pkin(4);
t48 = mrSges(3,2) - mrSges(4,3);
t46 = -m(4) - m(5) - m(6);
t45 = t42 - t49;
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t25 = t18 * mrSges(4,1) - t15 * mrSges(4,2);
t43 = t53 * t15 + mrSges(3,1) + t25;
t16 = sin(qJ(1));
t38 = t16 * pkin(1);
t19 = cos(qJ(1));
t11 = t19 * pkin(1);
t37 = t9 * t14;
t36 = t10 * t14;
t35 = t14 * t18;
t32 = t18 * t17;
t27 = t18 * pkin(3) + t15 * pkin(7);
t26 = -t15 * t13 + t18 * t8;
t1 = t10 * t17 + t9 * t35;
t3 = -t10 * t35 + t9 * t17;
t4 = t10 * t32 + t37;
t2 = -t9 * t32 + t36;
t5 = [(-t37 * t42 - m(3) * t11 - t19 * mrSges(2,1) + t16 * mrSges(2,2) + t48 * t9 + t49 * t4 + t46 * (t10 * pkin(2) + t9 * pkin(6) + t11) - t47 * t3 + (-m(5) * t27 - m(6) * t26 - t43) * t10) * g(2) + (-t36 * t42 + m(3) * t38 + t16 * mrSges(2,1) + t19 * mrSges(2,2) + t46 * (t10 * pkin(6) - t38) + t49 * t2 + t48 * t10 - t47 * t1 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t27) - m(6) * (-pkin(2) - t26) + t43) * t9) * g(1), (-m(3) + t46) * g(3), -g(3) * t25 + (t52 * g(3) + t50 * (mrSges(4,2) + t51)) * t18 + (t51 * g(3) + t50 * (mrSges(4,1) - t52)) * t15, (t45 * t14 + t47 * t17) * g(3) * t15 + (t45 * t1 - t47 * t2) * g(2) + (-t45 * t3 + t47 * t4) * g(1), (g(3) * t18 - t50 * t15) * m(6)];
taug = t5(:);
