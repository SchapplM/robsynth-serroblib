% Calculate Gravitation load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:56
% EndTime: 2020-01-03 11:35:58
% DurationCPUTime: 0.27s
% Computational Cost: add. (261->60), mult. (185->76), div. (0->0), fcn. (158->10), ass. (0->37)
t34 = sin(pkin(9));
t61 = t34 * mrSges(5,2) - mrSges(4,1);
t60 = mrSges(4,2) - mrSges(5,3);
t59 = m(5) + m(6);
t33 = qJ(1) + pkin(8);
t30 = qJ(3) + t33;
t26 = sin(t30);
t35 = cos(pkin(9));
t55 = t26 * t35;
t56 = t26 * t34;
t58 = pkin(4) * t55 + pkin(7) * t56;
t27 = cos(t30);
t53 = t27 * t35;
t54 = t27 * t34;
t57 = pkin(4) * t53 + pkin(7) * t54;
t36 = sin(qJ(5));
t51 = t35 * t36;
t38 = cos(qJ(5));
t50 = t35 * t38;
t49 = t27 * pkin(3) + t26 * qJ(4);
t28 = sin(t33);
t37 = sin(qJ(1));
t48 = t37 * pkin(1) + pkin(2) * t28;
t29 = cos(t33);
t39 = cos(qJ(1));
t47 = t39 * pkin(1) + pkin(2) * t29;
t46 = -m(3) * pkin(1) - mrSges(2,1);
t45 = t26 * pkin(3) - qJ(4) * t27;
t44 = t47 + t49;
t42 = t45 + t48;
t5 = -t26 * t51 - t27 * t38;
t6 = t26 * t50 - t27 * t36;
t41 = -mrSges(5,1) * t55 - t6 * mrSges(6,1) - t5 * mrSges(6,2) - mrSges(6,3) * t56 + t61 * t26 - t60 * t27;
t7 = -t26 * t38 + t27 * t51;
t8 = t26 * t36 + t27 * t50;
t40 = -mrSges(5,1) * t53 - t8 * mrSges(6,1) + t7 * mrSges(6,2) - mrSges(6,3) * t54 + t60 * t26 + t61 * t27;
t1 = [(-mrSges(2,2) * t39 - mrSges(3,1) * t28 - mrSges(3,2) * t29 - m(4) * t48 - m(5) * t42 - m(6) * (t42 + t58) + t46 * t37 + t41) * g(3) + (t37 * mrSges(2,2) - t29 * mrSges(3,1) + t28 * mrSges(3,2) - m(4) * t47 - m(5) * t44 - m(6) * (t44 + t57) + t46 * t39 + t40) * g(2), (-m(3) - m(4) - t59) * g(1), (-m(5) * t45 - m(6) * (t45 + t58) + t41) * g(3) + (-m(5) * t49 - m(6) * (t49 + t57) + t40) * g(2), t59 * (g(2) * t27 + g(3) * t26), -g(2) * (mrSges(6,1) * t5 - mrSges(6,2) * t6) - g(3) * (mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(1) * (-mrSges(6,1) * t36 - mrSges(6,2) * t38) * t34];
taug = t1(:);
