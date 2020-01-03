% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:04
% DurationCPUTime: 0.26s
% Computational Cost: add. (348->64), mult. (223->67), div. (0->0), fcn. (169->10), ass. (0->42)
t34 = qJ(4) + qJ(5);
t27 = sin(t34);
t29 = cos(t34);
t68 = t29 * mrSges(6,1) - t27 * mrSges(6,2);
t36 = sin(qJ(4));
t67 = t36 * mrSges(5,2) - t68;
t66 = mrSges(6,1) * t27 + mrSges(6,2) * t29;
t65 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t38 = cos(qJ(4));
t64 = -t38 * mrSges(5,1) - mrSges(4,1) + t67;
t49 = m(6) * pkin(4) + mrSges(5,1);
t63 = -mrSges(5,2) * t38 - t49 * t36;
t35 = qJ(1) + qJ(2);
t31 = qJ(3) + t35;
t24 = sin(t31);
t25 = cos(t31);
t26 = t38 * pkin(4) + pkin(3);
t40 = -pkin(9) - pkin(8);
t62 = t24 * t26 + t25 * t40;
t61 = t66 * t25;
t28 = sin(t35);
t22 = pkin(2) * t28;
t30 = cos(t35);
t23 = pkin(2) * t30;
t60 = g(2) * t24;
t54 = t25 * pkin(3) + t24 * pkin(8);
t37 = sin(qJ(1));
t32 = t37 * pkin(1);
t53 = t22 + t32;
t52 = t22 + t62;
t51 = t23 + t54;
t50 = -m(3) * pkin(1) - mrSges(2,1);
t48 = -t24 * t40 + t25 * t26;
t47 = t23 + t48;
t44 = t65 * t24 + t64 * t25;
t43 = (m(5) * pkin(8) - t65) * t25 + t64 * t24;
t42 = -t30 * mrSges(3,1) + t28 * mrSges(3,2) + t44;
t41 = -t28 * mrSges(3,1) - t30 * mrSges(3,2) + t43;
t39 = cos(qJ(1));
t33 = t39 * pkin(1);
t16 = t24 * pkin(3);
t1 = [(-t39 * mrSges(2,2) - m(4) * t53 - m(5) * (t16 + t53) - m(6) * (t32 + t52) + t50 * t37 + t41) * g(3) + (t37 * mrSges(2,2) - m(4) * (t23 + t33) - m(5) * (t33 + t51) - m(6) * (t33 + t47) + t50 * t39 + t42) * g(2), (-m(4) * t22 - m(5) * (t16 + t22) - m(6) * t52 + t41) * g(3) + (-m(4) * t23 - m(5) * t51 - m(6) * t47 + t42) * g(2), (-m(5) * t16 - m(6) * t62 + t43) * g(3) + (-m(5) * t54 - m(6) * t48 + t44) * g(2), (t63 * t25 - t61) * g(3) + (-t49 * t38 + t67) * g(1) + (t66 - t63) * t60, -g(1) * t68 - g(3) * t61 + t60 * t66];
taug = t1(:);
