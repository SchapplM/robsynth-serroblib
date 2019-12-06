% Calculate Gravitation load on the joints for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:19
% EndTime: 2019-12-05 16:25:22
% DurationCPUTime: 0.73s
% Computational Cost: add. (316->76), mult. (591->117), div. (0->0), fcn. (662->12), ass. (0->38)
t28 = sin(qJ(5));
t31 = cos(qJ(5));
t68 = -m(6) * pkin(4) - mrSges(6,1) * t31 + mrSges(6,2) * t28 - mrSges(5,1);
t63 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t70 = m(5) + m(6);
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t53 = cos(pkin(5));
t26 = sin(pkin(5));
t30 = sin(qJ(2));
t56 = t26 * t30;
t73 = -t29 * t56 + t53 * t32;
t33 = cos(qJ(2));
t25 = sin(pkin(9));
t47 = t25 * t53;
t52 = cos(pkin(9));
t14 = -t30 * t47 + t52 * t33;
t55 = t26 * t32;
t72 = -t14 * t29 + t25 * t55;
t64 = -m(4) * pkin(7) - t28 * mrSges(6,1) - t31 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t24 = qJ(3) + pkin(10);
t22 = sin(t24);
t23 = cos(t24);
t71 = -m(4) * pkin(2) - t32 * mrSges(4,1) + t29 * mrSges(4,2) + t63 * t22 + t68 * t23 - mrSges(3,1);
t39 = t53 * t52;
t12 = t25 * t33 + t30 * t39;
t46 = t26 * t52;
t36 = -t12 * t29 - t32 * t46;
t57 = t25 * t26;
t54 = t26 * t33;
t27 = -qJ(4) - pkin(7);
t21 = t32 * pkin(3) + pkin(2);
t13 = t52 * t30 + t33 * t47;
t11 = t25 * t30 - t33 * t39;
t8 = t53 * t22 + t23 * t56;
t4 = t14 * t23 + t22 * t57;
t2 = t12 * t23 - t22 * t46;
t1 = [(-m(2) - m(3) - m(4) - t70) * g(3), (-t70 * (-t11 * t21 - t12 * t27) + t64 * t12 - t71 * t11) * g(2) + (-t70 * (-t13 * t21 - t14 * t27) + t64 * t14 - t71 * t13) * g(1) + (-t70 * t21 * t54 + (t71 * t33 + (t70 * t27 + t64) * t30) * t26) * g(3), (-t73 * mrSges(4,1) - (-t53 * t29 - t30 * t55) * mrSges(4,2) + t63 * t8 + t68 * (-t22 * t56 + t53 * t23)) * g(3) + (-t36 * mrSges(4,1) - (-t12 * t32 + t29 * t46) * mrSges(4,2) + t68 * (-t12 * t22 - t23 * t46) + t63 * t2) * g(2) + (-t72 * mrSges(4,1) - (-t14 * t32 - t29 * t57) * mrSges(4,2) + t63 * t4 + t68 * (-t14 * t22 + t23 * t57)) * g(1) + (-g(1) * t72 - g(2) * t36 - g(3) * t73) * t70 * pkin(3), t70 * (-g(1) * t13 - g(2) * t11 + g(3) * t54), -g(1) * ((t13 * t31 - t28 * t4) * mrSges(6,1) + (-t13 * t28 - t31 * t4) * mrSges(6,2)) - g(2) * ((t11 * t31 - t2 * t28) * mrSges(6,1) + (-t11 * t28 - t2 * t31) * mrSges(6,2)) - g(3) * ((-t8 * t28 - t31 * t54) * mrSges(6,1) + (t28 * t54 - t8 * t31) * mrSges(6,2))];
taug = t1(:);
