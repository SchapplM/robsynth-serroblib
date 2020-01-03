% Calculate Gravitation load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:05
% DurationCPUTime: 0.76s
% Computational Cost: add. (207->89), mult. (486->105), div. (0->0), fcn. (483->8), ass. (0->46)
t83 = mrSges(4,1) + mrSges(5,1);
t73 = mrSges(4,2) - mrSges(5,3);
t82 = mrSges(5,2) + mrSges(4,3);
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t72 = g(1) * t29 + g(2) * t26;
t81 = m(6) * pkin(7) + mrSges(6,3) - t82;
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t24 = sin(qJ(5));
t27 = cos(qJ(5));
t37 = t22 * t24 + t23 * t27;
t38 = t22 * t27 - t23 * t24;
t80 = t37 * mrSges(6,1) + t38 * mrSges(6,2) - t73 * t22 + t83 * t23;
t77 = pkin(3) * t23 + qJ(4) * t22;
t75 = m(5) + m(6);
t74 = mrSges(2,2) - mrSges(3,3);
t71 = m(4) + t75;
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t42 = t28 * mrSges(3,1) - t25 * mrSges(3,2);
t69 = t82 * t25 + t42;
t66 = m(6) * pkin(4) + t83;
t63 = pkin(4) * t23;
t60 = g(3) * t25;
t19 = t28 * pkin(2);
t56 = t26 * t28;
t55 = t28 * t29;
t54 = t29 * t22;
t17 = t25 * qJ(3);
t53 = t19 + t17;
t52 = t29 * pkin(1) + t26 * pkin(6);
t48 = -pkin(1) - t19;
t47 = t77 * t28 + t53;
t46 = pkin(2) * t55 + t29 * t17 + t52;
t5 = t22 * t56 + t23 * t29;
t6 = t23 * t56 - t54;
t44 = t24 * t6 - t27 * t5;
t43 = -t24 * t5 - t27 * t6;
t35 = -pkin(2) - t77;
t20 = t29 * pkin(6);
t8 = t26 * t22 + t23 * t55;
t7 = -t26 * t23 + t28 * t54;
t2 = t24 * t7 + t27 * t8;
t1 = -t24 * t8 + t27 * t7;
t3 = [(-m(3) * t52 - m(4) * t46 - t2 * mrSges(6,1) - t1 * mrSges(6,2) - t75 * (t8 * pkin(3) + t7 * qJ(4) + t46) - t66 * t8 + t73 * t7 + t74 * t26 + (t81 * t25 - mrSges(2,1) - t42) * t29) * g(2) + (-t43 * mrSges(6,1) - t44 * mrSges(6,2) - t75 * (-t6 * pkin(3) - qJ(4) * t5 + t20) + t66 * t6 - t73 * t5 + t74 * t29 + (-m(3) - m(4)) * t20 + (mrSges(2,1) + m(3) * pkin(1) - m(6) * t48 - (m(6) * (pkin(7) - qJ(3)) + mrSges(6,3)) * t25 + (-m(4) - m(5)) * (t48 - t17) + t69) * t26) * g(1), (-m(4) * t53 - m(5) * t47 - m(6) * (-pkin(7) * t25 + t47) + t25 * mrSges(6,3) + (-m(6) * t63 - t80) * t28 - t69) * g(3) + ((mrSges(3,1) - m(6) * (t35 - t63) - m(5) * t35 + m(4) * pkin(2) + t80) * t25 + (-qJ(3) * t71 + mrSges(3,2) + t81) * t28) * t72, (g(3) * t28 - t72 * t25) * t71, t75 * (-g(1) * t7 - g(2) * t5 - t22 * t60), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (-t44 * mrSges(6,1) + t43 * mrSges(6,2)) - (t38 * mrSges(6,1) - t37 * mrSges(6,2)) * t60];
taug = t3(:);
