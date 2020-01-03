% Calculate Gravitation load on the joints for
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:50
% EndTime: 2020-01-03 12:08:51
% DurationCPUTime: 0.25s
% Computational Cost: add. (289->69), mult. (227->74), div. (0->0), fcn. (173->10), ass. (0->43)
t40 = qJ(3) + pkin(9);
t31 = sin(t40);
t32 = cos(t40);
t43 = sin(qJ(3));
t33 = qJ(5) + t40;
t25 = sin(t33);
t64 = t25 * mrSges(6,2);
t71 = -t32 * mrSges(5,1) + t43 * mrSges(4,2) + t31 * mrSges(5,2) + t64;
t45 = cos(qJ(3));
t70 = -t45 * mrSges(4,1) - mrSges(3,1) + t71;
t69 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t55 = m(5) * pkin(3) + mrSges(4,1);
t68 = -t55 * t43 + m(6) * (-t43 * pkin(3) - pkin(4) * t31) - mrSges(5,1) * t31 - mrSges(4,2) * t45 - mrSges(5,2) * t32;
t41 = qJ(1) + qJ(2);
t34 = sin(t41);
t67 = g(2) * t34;
t35 = cos(t41);
t42 = -qJ(4) - pkin(7);
t39 = -pkin(8) + t42;
t37 = t45 * pkin(3);
t59 = pkin(4) * t32 + t37;
t7 = pkin(2) + t59;
t66 = t34 * t7 + t35 * t39;
t65 = mrSges(6,1) * t25;
t26 = cos(t33);
t17 = t26 * mrSges(6,1);
t63 = t26 * t35;
t30 = t37 + pkin(2);
t60 = t34 * t30 + t35 * t42;
t58 = t35 * pkin(2) + t34 * pkin(7);
t57 = g(3) * (mrSges(6,2) * t63 + t35 * t65);
t56 = -m(3) * pkin(1) - mrSges(2,1);
t54 = -t34 * t39 + t35 * t7;
t53 = t35 * t30 - t34 * t42;
t51 = -mrSges(6,2) * t26 - t65;
t48 = -mrSges(6,1) * t63 + t69 * t34 + t70 * t35;
t47 = (m(4) * pkin(7) - t69) * t35 + (-t17 + t70) * t34;
t46 = cos(qJ(1));
t44 = sin(qJ(1));
t38 = t46 * pkin(1);
t36 = t44 * pkin(1);
t28 = t34 * pkin(2);
t1 = [(-t46 * mrSges(2,2) - m(4) * (t28 + t36) - m(5) * (t36 + t60) - m(6) * (t36 + t66) + t56 * t44 + t47) * g(3) + (t44 * mrSges(2,2) - m(4) * (t38 + t58) - m(5) * (t38 + t53) - m(6) * (t38 + t54) + t56 * t46 + t48) * g(2), (-m(4) * t28 - m(5) * t60 - m(6) * t66 + t47) * g(3) + (-m(4) * t58 - m(5) * t53 - m(6) * t54 + t48) * g(2), -t57 + (-m(6) * t59 - t55 * t45 - t17 + t71) * g(1) + t68 * g(3) * t35 + (-t51 - t68) * t67, (m(5) + m(6)) * (g(2) * t35 + g(3) * t34), -g(1) * (t17 - t64) - t57 - t51 * t67];
taug = t1(:);
