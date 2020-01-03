% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:26
% EndTime: 2020-01-03 11:30:28
% DurationCPUTime: 0.37s
% Computational Cost: add. (216->72), mult. (272->85), div. (0->0), fcn. (249->10), ass. (0->44)
t62 = -m(5) - m(6);
t42 = m(4) - t62;
t61 = m(3) + t42;
t60 = -m(6) * pkin(4) - mrSges(5,1);
t29 = cos(pkin(9));
t17 = t29 * pkin(3) + pkin(2);
t26 = pkin(9) + qJ(4);
t19 = cos(t26);
t27 = sin(pkin(9));
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t31 = -pkin(6) - qJ(3);
t59 = -mrSges(2,1) + (-mrSges(3,1) - m(4) * pkin(2) - t29 * mrSges(4,1) + t27 * mrSges(4,2) - m(5) * t17 - m(6) * (pkin(4) * t19 + t17)) * t30 + (mrSges(3,2) - m(4) * qJ(3) - mrSges(4,3) + m(5) * t31 - mrSges(5,3) - m(6) * (pkin(7) - t31) - mrSges(6,3)) * t28;
t52 = t27 * pkin(3);
t18 = sin(t26);
t54 = pkin(4) * t18;
t58 = -m(5) * t52 - m(6) * (t52 + t54) + mrSges(2,2) - mrSges(3,3) - t27 * mrSges(4,1) - t29 * mrSges(4,2);
t20 = qJ(5) + t26;
t16 = cos(t20);
t33 = cos(qJ(1));
t46 = t33 * t16;
t15 = sin(t20);
t32 = sin(qJ(1));
t51 = t32 * t15;
t5 = -t30 * t51 - t46;
t47 = t33 * t15;
t50 = t32 * t16;
t6 = t30 * t50 - t47;
t56 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t7 = t30 * t47 - t50;
t8 = t30 * t46 + t51;
t55 = t7 * mrSges(6,1) + t8 * mrSges(6,2);
t53 = g(1) * t28;
t49 = t32 * t18;
t48 = t32 * t19;
t45 = t33 * t18;
t44 = t33 * t19;
t37 = -mrSges(6,1) * t15 - mrSges(6,2) * t16;
t11 = t30 * t45 - t48;
t9 = -t30 * t49 - t44;
t23 = t32 * pkin(1);
t12 = t30 * t44 + t49;
t10 = t30 * t48 - t45;
t1 = [(-m(4) * t23 - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (-m(3) + t62) * (-t33 * qJ(2) + t23) + (m(4) * qJ(2) - t58) * t33 + t59 * t32) * g(3) + (-t12 * mrSges(5,1) - t8 * mrSges(6,1) + t11 * mrSges(5,2) + t7 * mrSges(6,2) - t61 * (t33 * pkin(1) + t32 * qJ(2)) + t58 * t32 + t59 * t33) * g(2), (g(2) * t33 + g(3) * t32) * t61, (t30 * g(1) + (-g(2) * t32 + g(3) * t33) * t28) * t42, (m(6) * t54 + mrSges(5,1) * t18 + mrSges(5,2) * t19 - t37) * t53 + (-mrSges(5,2) * t12 + t60 * t11 - t55) * g(3) + (mrSges(5,2) * t10 + t60 * t9 - t56) * g(2), -g(2) * t56 - g(3) * t55 - t37 * t53];
taug = t1(:);
