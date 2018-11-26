% Calculate Gravitation load on the joints for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:56:22
% EndTime: 2018-11-23 15:56:22
% DurationCPUTime: 0.59s
% Computational Cost: add. (384->79), mult. (351->87), div. (0->0), fcn. (301->10), ass. (0->45)
t58 = -mrSges(6,1) - mrSges(7,1);
t57 = mrSges(6,2) + mrSges(7,2);
t15 = qJ(1) + pkin(9);
t11 = cos(t15);
t9 = sin(t15);
t60 = g(1) * t11 + g(2) * t9;
t14 = qJ(3) + pkin(10);
t8 = sin(t14);
t64 = t60 * t8;
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t6 = t21 * pkin(5) + pkin(4);
t63 = -m(6) * pkin(4) - m(7) * t6 + t57 * t18 + t58 * t21;
t37 = m(5) + m(6) + m(7);
t62 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t10 = cos(t14);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t61 = -t22 * mrSges(4,1) - t10 * mrSges(5,1) + t19 * mrSges(4,2) - t62 * t8;
t51 = m(7) * pkin(5);
t59 = m(3) + m(4);
t56 = t51 - t58;
t53 = m(4) * pkin(2) + mrSges(3,1) - t61;
t52 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t48 = t8 * pkin(8);
t20 = sin(qJ(1));
t45 = t20 * pkin(1);
t12 = t22 * pkin(3);
t23 = cos(qJ(1));
t13 = t23 * pkin(1);
t16 = -qJ(6) - pkin(8);
t42 = t8 * t16;
t41 = t9 * t18;
t40 = t9 * t21;
t39 = t11 * t18;
t38 = t11 * t21;
t35 = t10 * pkin(4) + t48;
t34 = t10 * t6 - t42;
t1 = t10 * t41 + t38;
t3 = -t10 * t39 + t40;
t17 = -qJ(4) - pkin(7);
t7 = t12 + pkin(2);
t4 = t10 * t38 + t41;
t2 = -t10 * t40 + t39;
t5 = [(-t41 * t51 - t23 * mrSges(2,1) + t20 * mrSges(2,2) + t58 * t4 - t37 * (t11 * t7 - t9 * t17 + t13) - t57 * t3 - t59 * t13 + t52 * t9 + (-m(6) * t35 - m(7) * t34 - t53) * t11) * g(2) + (-t39 * t51 + t20 * mrSges(2,1) + t23 * mrSges(2,2) + t59 * t45 - t37 * (-t11 * t17 - t45) + t58 * t2 - t57 * t1 + t52 * t11 + (m(5) * t7 - m(6) * (-t35 - t7) - m(7) * (-t34 - t7) + t53) * t9) * g(1) (-t37 - t59) * g(3) (mrSges(5,1) - t63) * t64 + (-m(5) * t12 - m(6) * (t12 + t48) - m(7) * (t12 - t42) + t61 + t63 * t10) * g(3) + (mrSges(4,2) * t22 + (-m(6) * pkin(8) + m(7) * t16 - t62) * t10 + (t37 * pkin(3) + mrSges(4,1)) * t19) * t60 (-g(1) * t9 + g(2) * t11) * t37 (t56 * t18 + t57 * t21) * g(3) * t8 + (t56 * t1 - t57 * t2) * g(2) + (-t56 * t3 + t57 * t4) * g(1) (g(3) * t10 - t64) * m(7)];
taug  = t5(:);
