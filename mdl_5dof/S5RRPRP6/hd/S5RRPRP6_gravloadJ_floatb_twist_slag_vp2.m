% Calculate Gravitation load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:56:48
% DurationCPUTime: 0.52s
% Computational Cost: add. (242->68), mult. (328->77), div. (0->0), fcn. (287->8), ass. (0->39)
t54 = -mrSges(5,1) - mrSges(6,1);
t53 = mrSges(5,2) + mrSges(6,2);
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t55 = g(1) * t19 + g(2) * t16;
t11 = qJ(2) + pkin(8);
t8 = sin(t11);
t59 = t55 * t8;
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t6 = t17 * pkin(4) + pkin(3);
t58 = -m(5) * pkin(3) - m(6) * t6 + t53 * t14 + t54 * t17;
t52 = m(4) + m(5) + m(6);
t57 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t15 = sin(qJ(2));
t18 = cos(qJ(2));
t9 = cos(t11);
t56 = -t18 * mrSges(3,1) - t9 * mrSges(4,1) + t15 * mrSges(3,2) - t57 * t8;
t46 = m(6) * pkin(4);
t51 = t46 - t54;
t50 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t47 = m(3) * pkin(1) + mrSges(2,1) - t56;
t44 = t8 * pkin(7);
t10 = t18 * pkin(2);
t12 = -qJ(5) - pkin(7);
t38 = t8 * t12;
t37 = t16 * t14;
t36 = t16 * t17;
t35 = t19 * t14;
t34 = t19 * t17;
t30 = t9 * pkin(3) + t44;
t29 = t9 * t6 - t38;
t3 = -t9 * t35 + t36;
t1 = t9 * t37 + t34;
t13 = -qJ(3) - pkin(6);
t7 = t10 + pkin(1);
t4 = t9 * t34 + t37;
t2 = -t9 * t36 + t35;
t5 = [(-t37 * t46 + t54 * t4 - t52 * (-t16 * t13 + t19 * t7) - t53 * t3 + t50 * t16 + (-m(5) * t30 - m(6) * t29 - t47) * t19) * g(2) + (t54 * t2 - t53 * t1 + (t52 * t13 - t46 * t14 + t50) * t19 + (m(4) * t7 - m(5) * (-t30 - t7) - m(6) * (-t29 - t7) + t47) * t16) * g(1), (mrSges(4,1) - t58) * t59 + (-m(4) * t10 - m(5) * (t10 + t44) - m(6) * (t10 - t38) + t56 + t58 * t9) * g(3) + (mrSges(3,2) * t18 + (-m(5) * pkin(7) + m(6) * t12 - t57) * t9 + (t52 * pkin(2) + mrSges(3,1)) * t15) * t55, (-g(1) * t16 + g(2) * t19) * t52, (t51 * t14 + t53 * t17) * g(3) * t8 + (t51 * t1 - t53 * t2) * g(2) + (-t51 * t3 + t53 * t4) * g(1), (g(3) * t9 - t59) * m(6)];
taug = t5(:);
