% Calculate Gravitation load on the joints for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:56
% EndTime: 2019-12-05 18:52:57
% DurationCPUTime: 0.24s
% Computational Cost: add. (271->64), mult. (272->80), div. (0->0), fcn. (230->10), ass. (0->44)
t38 = cos(qJ(3));
t61 = m(5) + m(6);
t48 = t61 * pkin(2);
t47 = mrSges(4,1) + t48;
t66 = t47 * t38;
t35 = sin(qJ(3));
t53 = t35 * mrSges(4,2);
t32 = qJ(3) + qJ(4);
t28 = sin(t32);
t59 = t28 * mrSges(5,2);
t46 = t53 + t59;
t65 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t21 = t28 * mrSges(6,3);
t30 = cos(t32);
t64 = t30 * mrSges(5,1) + t21;
t63 = mrSges(2,1) + (m(3) + m(4)) * pkin(1);
t60 = pkin(2) * t38;
t33 = qJ(1) + qJ(2);
t31 = cos(t33);
t58 = t30 * t31;
t34 = sin(qJ(5));
t57 = t30 * t34;
t37 = cos(qJ(5));
t56 = t30 * t37;
t55 = t31 * t34;
t54 = t31 * t37;
t52 = t38 * mrSges(4,1);
t29 = sin(t33);
t49 = mrSges(6,2) * t28 * t34;
t51 = (mrSges(6,3) * t30 + t49) * t29;
t50 = mrSges(6,3) * t58 + t31 * t49;
t45 = mrSges(3,1) + t64;
t5 = t29 * t57 + t54;
t6 = -t29 * t56 + t55;
t44 = -t6 * mrSges(6,1) - t5 * mrSges(6,2) - t46 * t29 + t65 * t31;
t43 = mrSges(5,2) * t30 + (mrSges(6,1) * t37 + mrSges(5,1)) * t28;
t42 = -mrSges(6,1) * t56 + mrSges(6,2) * t57 + t59 - t64;
t7 = t29 * t37 - t30 * t55;
t8 = t29 * t34 + t30 * t54;
t41 = -mrSges(5,1) * t58 - t8 * mrSges(6,1) - t7 * mrSges(6,2) + (-mrSges(3,1) - t21 - t52) * t31 + t65 * t29;
t40 = mrSges(4,2) * t38 + t47 * t35 + t43;
t39 = cos(qJ(1));
t36 = sin(qJ(1));
t1 = [(t36 * mrSges(2,2) + t46 * t31 - t61 * (t39 * pkin(1) + t31 * t60) - t63 * t39 + t41) * g(2) + (t39 * mrSges(2,2) - t61 * (-t36 * pkin(1) - t29 * t60) + t63 * t36 + (t45 + t52) * t29 + t44) * g(1), ((-t38 * t48 + t46) * t31 + t41) * g(2) + ((t45 + t66) * t29 + t44) * g(1), (t42 + t53 - t66) * g(3) + (t40 * t29 - t51) * g(2) + (t40 * t31 - t50) * g(1), t42 * g(3) + (t43 * t29 - t51) * g(2) + (t43 * t31 - t50) * g(1), -g(1) * (mrSges(6,1) * t7 - mrSges(6,2) * t8) - g(2) * (-t5 * mrSges(6,1) + mrSges(6,2) * t6) - g(3) * (-mrSges(6,1) * t34 - mrSges(6,2) * t37) * t28];
taug = t1(:);
