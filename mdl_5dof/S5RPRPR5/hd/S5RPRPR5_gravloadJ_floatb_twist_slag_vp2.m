% Calculate Gravitation load on the joints for
% S5RPRPR5
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:52
% EndTime: 2019-12-05 17:55:54
% DurationCPUTime: 0.41s
% Computational Cost: add. (233->87), mult. (307->105), div. (0->0), fcn. (285->10), ass. (0->54)
t62 = m(5) + m(6);
t66 = m(3) + m(4) + t62;
t70 = -t66 * qJ(2) + mrSges(2,2) - mrSges(3,3);
t63 = m(5) * pkin(3);
t68 = mrSges(4,1) + t63;
t30 = qJ(3) + pkin(9);
t25 = cos(t30);
t36 = cos(qJ(3));
t28 = t36 * pkin(3);
t19 = pkin(4) * t25 + t28;
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t33 = -qJ(4) - pkin(6);
t65 = mrSges(2,1) + m(3) * pkin(1) + t32 * mrSges(3,1) - m(4) * (-pkin(2) * t32 - pkin(1)) - m(5) * (-(t28 + pkin(2)) * t32 - pkin(1)) - m(6) * (-(pkin(2) + t19) * t32 - pkin(1)) + (-mrSges(3,2) + m(4) * pkin(6) + mrSges(4,3) - m(5) * t33 + mrSges(5,3) - m(6) * (-pkin(7) + t33) + mrSges(6,3)) * t31;
t24 = sin(t30);
t34 = sin(qJ(3));
t58 = t34 * pkin(3);
t18 = pkin(4) * t24 + t58;
t64 = m(5) * t58 + m(6) * t18;
t26 = qJ(5) + t30;
t22 = cos(t26);
t37 = cos(qJ(1));
t49 = t37 * t22;
t21 = sin(t26);
t35 = sin(qJ(1));
t57 = t35 * t21;
t5 = t32 * t57 + t49;
t50 = t37 * t21;
t56 = t35 * t22;
t6 = t32 * t56 - t50;
t61 = t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = t32 * t50 - t56;
t8 = -t32 * t49 - t57;
t60 = -t7 * mrSges(6,1) + t8 * mrSges(6,2);
t59 = g(1) * t31;
t55 = t35 * t24;
t54 = t35 * t25;
t53 = t35 * t34;
t52 = t35 * t36;
t51 = t37 * t18;
t48 = t37 * t24;
t47 = t37 * t25;
t46 = t37 * t34;
t45 = t37 * t36;
t42 = -mrSges(6,1) * t21 - mrSges(6,2) * t22;
t15 = t32 * t46 - t52;
t13 = t32 * t53 + t45;
t16 = -t32 * t45 - t53;
t14 = t32 * t52 - t46;
t12 = -t32 * t47 - t55;
t11 = t32 * t48 - t54;
t10 = t32 * t54 - t48;
t9 = t32 * t55 + t47;
t1 = [(-m(6) * t51 + t14 * mrSges(4,1) + t10 * mrSges(5,1) + t6 * mrSges(6,1) - t13 * mrSges(4,2) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + t65 * t35 + t70 * t37 - t46 * t63) * g(3) + (-t16 * mrSges(4,1) - t12 * mrSges(5,1) - t8 * mrSges(6,1) - t15 * mrSges(4,2) - t11 * mrSges(5,2) - t7 * mrSges(6,2) + t65 * t37 + (t64 - t70) * t35) * g(2), -(g(2) * t37 + g(3) * t35) * t66, (mrSges(4,1) * t34 + mrSges(5,1) * t24 + mrSges(4,2) * t36 + mrSges(5,2) * t25 - t42 + t64) * t59 + (-t16 * mrSges(4,2) + t11 * mrSges(5,1) - t12 * mrSges(5,2) - m(6) * (t35 * t19 - t32 * t51) - t60 + t68 * t15) * g(3) + (-t14 * mrSges(4,2) - t9 * mrSges(5,1) - t10 * mrSges(5,2) - m(6) * (t35 * t32 * t18 + t37 * t19) - t61 - t68 * t13) * g(2), (t32 * g(1) + t31 * (g(2) * t35 - g(3) * t37)) * t62, -g(2) * t61 - g(3) * t60 - t42 * t59];
taug = t1(:);
