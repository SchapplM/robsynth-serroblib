% Calculate Gravitation load on the joints for
% S5RRPPR7
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:34:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (212->73), mult. (279->78), div. (0->0), fcn. (233->8), ass. (0->40)
t51 = m(5) + m(6);
t65 = m(4) + t51;
t64 = mrSges(4,1) - mrSges(5,2);
t63 = -mrSges(4,2) + mrSges(5,3);
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t55 = g(1) * t23 + g(2) * t20;
t16 = qJ(2) + pkin(8);
t13 = sin(t16);
t14 = cos(t16);
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t60 = -t22 * mrSges(3,1) + t19 * mrSges(3,2) - t63 * t13 - t64 * t14;
t10 = t13 * qJ(4);
t43 = t23 * t14;
t58 = pkin(3) * t43 + t23 * t10;
t53 = -m(3) * pkin(1) - mrSges(2,1) + t60;
t17 = -qJ(3) - pkin(6);
t52 = -m(3) * pkin(6) - m(6) * (pkin(4) - t17) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t47 = g(3) * t14;
t11 = t14 * pkin(3);
t15 = t22 * pkin(2);
t18 = sin(qJ(5));
t45 = t20 * t18;
t21 = cos(qJ(5));
t44 = t20 * t21;
t42 = t23 * t18;
t41 = t23 * t21;
t38 = t11 + t10 + t15;
t12 = t15 + pkin(1);
t9 = t23 * t12;
t36 = -t20 * t17 + t9;
t35 = -t12 - t10;
t34 = m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3);
t29 = t18 * mrSges(6,1) + t21 * mrSges(6,2);
t4 = -t13 * t45 + t41;
t3 = t13 * t44 + t42;
t2 = t13 * t42 + t44;
t1 = t13 * t41 - t45;
t5 = [(-m(4) * t36 - m(5) * (t36 + t58) - m(6) * (pkin(7) * t43 + t58 + t9) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - mrSges(6,3) * t43 + t53 * t23 + t52 * t20) * g(2) + (-t4 * mrSges(6,1) + t3 * mrSges(6,2) + ((m(4) + m(5)) * t17 + t52) * t23 + (m(4) * t12 - m(5) * (t35 - t11) - m(6) * t35 - t14 * t34 - t53) * t20) * g(1), (-m(4) * t15 - m(5) * t38 - m(6) * (t14 * pkin(7) + t38) - t14 * mrSges(6,3) - t29 * t13 + t60) * g(3) + (mrSges(3,2) * t22 + (m(5) * pkin(3) - t34 + t64) * t13 + (-qJ(4) * t51 - t29 - t63) * t14 + (t65 * pkin(2) + mrSges(3,1)) * t19) * t55, (-g(1) * t20 + g(2) * t23) * t65, (-t13 * t55 + t47) * t51, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - (-mrSges(6,1) * t21 + mrSges(6,2) * t18) * t47];
taug = t5(:);
