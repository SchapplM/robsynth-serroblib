% Calculate Gravitation load on the joints for
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:04
% DurationCPUTime: 0.51s
% Computational Cost: add. (286->85), mult. (343->97), div. (0->0), fcn. (306->10), ass. (0->51)
t27 = cos(qJ(4));
t14 = t27 * pkin(4) + pkin(3);
t22 = qJ(4) + qJ(5);
t18 = sin(t22);
t19 = cos(t22);
t24 = sin(qJ(4));
t75 = -m(5) * pkin(3) - m(6) * t14 - t27 * mrSges(5,1) - t19 * mrSges(6,1) + t24 * mrSges(5,2) + t18 * mrSges(6,2);
t70 = m(4) + m(5) + m(6);
t74 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t65 = m(6) * pkin(4);
t21 = qJ(2) + pkin(9);
t16 = sin(t21);
t17 = cos(t21);
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t73 = -t28 * mrSges(3,1) - t17 * mrSges(4,1) + t25 * mrSges(3,2) - t74 * t16;
t72 = t24 * t65;
t71 = mrSges(5,1) + t65;
t69 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t66 = m(3) * pkin(1) + mrSges(2,1) - t73;
t29 = cos(qJ(1));
t48 = t29 * t19;
t26 = sin(qJ(1));
t53 = t26 * t18;
t5 = t17 * t53 + t48;
t49 = t29 * t18;
t52 = t26 * t19;
t6 = -t17 * t52 + t49;
t64 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t17 * t49 + t52;
t8 = t17 * t48 + t53;
t63 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t58 = g(3) * t16;
t57 = t16 * pkin(7);
t20 = t28 * pkin(2);
t30 = -pkin(8) - pkin(7);
t54 = t16 * t30;
t51 = t26 * t24;
t50 = t26 * t27;
t47 = t29 * t24;
t46 = t29 * t27;
t42 = t17 * pkin(3) + t57;
t38 = -mrSges(6,1) * t18 - mrSges(6,2) * t19;
t37 = t17 * t14 - t54;
t11 = -t17 * t47 + t50;
t9 = t17 * t51 + t46;
t23 = -qJ(3) - pkin(6);
t15 = t20 + pkin(1);
t12 = t17 * t46 + t51;
t10 = -t17 * t50 + t47;
t1 = [(-t51 * t65 - t12 * mrSges(5,1) - t8 * mrSges(6,1) - t11 * mrSges(5,2) - t7 * mrSges(6,2) - t70 * (t29 * t15 - t26 * t23) + t69 * t26 + (-m(5) * t42 - m(6) * t37 - t66) * t29) * g(2) + (-t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (t70 * t23 + t69 - t72) * t29 + (m(4) * t15 - m(5) * (-t15 - t42) - m(6) * (-t15 - t37) + t66) * t26) * g(1), (-m(4) * t20 - m(5) * (t20 + t57) - m(6) * (t20 - t54) + t75 * t17 + t73) * g(3) + (g(1) * t29 + g(2) * t26) * (mrSges(3,2) * t28 + (-m(5) * pkin(7) + m(6) * t30 - t74) * t17 + (mrSges(4,1) - t75) * t16 + (t70 * pkin(2) + mrSges(3,1)) * t25), (-g(1) * t26 + g(2) * t29) * t70, (mrSges(5,1) * t24 + mrSges(5,2) * t27 - t38 + t72) * t58 + (-t10 * mrSges(5,2) + t71 * t9 - t64) * g(2) + (t12 * mrSges(5,2) - t71 * t11 - t63) * g(1), -g(1) * t63 - g(2) * t64 - t38 * t58];
taug = t1(:);
