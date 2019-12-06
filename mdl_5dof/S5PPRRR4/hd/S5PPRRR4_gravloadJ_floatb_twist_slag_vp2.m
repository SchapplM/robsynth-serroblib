% Calculate Gravitation load on the joints for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:36
% EndTime: 2019-12-05 15:18:38
% DurationCPUTime: 0.45s
% Computational Cost: add. (435->60), mult. (1210->103), div. (0->0), fcn. (1520->14), ass. (0->48)
t71 = m(5) + m(6);
t27 = sin(qJ(5));
t30 = cos(qJ(5));
t73 = m(6) * pkin(4) + t30 * mrSges(6,1) - t27 * mrSges(6,2) + mrSges(5,1);
t69 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t72 = pkin(3) * t71 - t69 * t28 + t73 * t31 + mrSges(4,1);
t53 = sin(pkin(11));
t54 = sin(pkin(10));
t42 = t54 * t53;
t57 = cos(pkin(11));
t58 = cos(pkin(10));
t49 = t58 * t57;
t60 = cos(pkin(5));
t35 = -t60 * t49 + t42;
t55 = sin(pkin(6));
t56 = sin(pkin(5));
t46 = t56 * t55;
t59 = cos(pkin(6));
t68 = t35 * t59 + t58 * t46;
t43 = t54 * t57;
t47 = t58 * t53;
t36 = t60 * t43 + t47;
t45 = t56 * t54;
t67 = t36 * t59 - t55 * t45;
t66 = t57 * t59 * t56 + t60 * t55;
t63 = m(3) + m(4) + t71;
t62 = -t27 * mrSges(6,1) - t30 * mrSges(6,2) - t71 * pkin(8) + mrSges(4,2) - mrSges(5,3);
t61 = cos(qJ(3));
t48 = t58 * t56;
t44 = t56 * t53;
t29 = sin(qJ(3));
t22 = -t60 * t42 + t49;
t21 = t60 * t47 + t43;
t20 = -t57 * t46 + t60 * t59;
t17 = t36 * t55 + t59 * t45;
t16 = t35 * t55 - t59 * t48;
t15 = t66 * t29 + t61 * t44;
t14 = t29 * t44 - t66 * t61;
t12 = t15 * t31 + t20 * t28;
t10 = t22 * t61 - t67 * t29;
t9 = t22 * t29 + t67 * t61;
t8 = t21 * t61 - t68 * t29;
t7 = t21 * t29 + t68 * t61;
t4 = t10 * t31 + t17 * t28;
t2 = t16 * t28 + t8 * t31;
t1 = [(-m(2) - t63) * g(3), (-g(1) * t45 + g(2) * t48 - g(3) * t60) * t63, (t72 * t14 + t62 * t15) * g(3) + (t62 * t8 + t72 * t7) * g(2) + (t62 * t10 + t72 * t9) * g(1), (t69 * t12 - t73 * (-t15 * t28 + t20 * t31)) * g(3) + (t69 * t2 - t73 * (t16 * t31 - t8 * t28)) * g(2) + (t69 * t4 - t73 * (-t10 * t28 + t17 * t31)) * g(1), -g(1) * ((-t4 * t27 + t9 * t30) * mrSges(6,1) + (-t9 * t27 - t4 * t30) * mrSges(6,2)) - g(2) * ((-t2 * t27 + t7 * t30) * mrSges(6,1) + (-t2 * t30 - t7 * t27) * mrSges(6,2)) - g(3) * ((-t12 * t27 + t14 * t30) * mrSges(6,1) + (-t12 * t30 - t14 * t27) * mrSges(6,2))];
taug = t1(:);
