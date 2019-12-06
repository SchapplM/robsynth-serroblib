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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:42:18
% EndTime: 2019-12-05 18:42:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (289->55), mult. (227->54), div. (0->0), fcn. (173->10), ass. (0->34)
t30 = qJ(3) + pkin(9);
t23 = sin(t30);
t24 = cos(t30);
t33 = sin(qJ(3));
t63 = t24 * mrSges(5,1) - t33 * mrSges(4,2) - t23 * mrSges(5,2);
t25 = qJ(5) + t30;
t20 = cos(t25);
t13 = t20 * mrSges(6,1);
t35 = cos(qJ(3));
t28 = t35 * pkin(3);
t47 = pkin(4) * t24 + t28;
t61 = m(4) * pkin(2) + m(5) * (t28 + pkin(2)) + m(6) * (pkin(2) + t47) + t35 * mrSges(4,1) + mrSges(3,1) + t13 + t63;
t32 = -qJ(4) - pkin(7);
t60 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-pkin(8) + t32) + m(5) * t32;
t58 = m(5) + m(6);
t45 = m(5) * pkin(3) + mrSges(4,1);
t55 = -t45 * t33 + m(6) * (-t33 * pkin(3) - pkin(4) * t23) - mrSges(5,1) * t23 - mrSges(4,2) * t35 - mrSges(5,2) * t24;
t31 = qJ(1) + qJ(2);
t27 = cos(t31);
t54 = g(3) * t27;
t53 = mrSges(6,2) * t20;
t19 = sin(t25);
t52 = t19 * mrSges(6,2);
t26 = sin(t31);
t51 = t19 * t26;
t46 = g(2) * (mrSges(6,1) * t51 + t26 * t53);
t44 = t13 - t52;
t42 = -mrSges(6,1) * t19 - t53;
t41 = mrSges(2,1) + (m(3) + m(4) + t58) * pkin(1);
t38 = (-t52 + t61) * t27 - t60 * t26;
t37 = -mrSges(6,2) * t51 + t61 * t26 + t60 * t27;
t36 = cos(qJ(1));
t34 = sin(qJ(1));
t1 = [(t36 * mrSges(2,2) + t41 * t34 + t37) * g(3) + (-t34 * mrSges(2,2) + t41 * t36 + t38) * g(2), t38 * g(2) + t37 * g(3), -t46 + (-m(6) * t47 - t45 * t35 - t44 - t63) * g(1) + t55 * g(2) * t26 + (-t42 - t55) * t54, t58 * (-g(2) * t27 - g(3) * t26), -g(1) * t44 - t42 * t54 - t46];
taug = t1(:);
