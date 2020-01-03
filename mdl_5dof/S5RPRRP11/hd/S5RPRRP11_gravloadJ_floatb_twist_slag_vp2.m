% Calculate Gravitation load on the joints for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:53:29
% DurationCPUTime: 0.50s
% Computational Cost: add. (251->66), mult. (337->78), div. (0->0), fcn. (310->8), ass. (0->39)
t62 = mrSges(5,1) + mrSges(6,1);
t61 = -mrSges(5,2) + mrSges(6,3);
t60 = -mrSges(5,3) - mrSges(6,2);
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t59 = t61 * t21 + t62 * t23;
t56 = m(5) + m(6);
t58 = m(4) + t56;
t19 = cos(pkin(8));
t57 = -mrSges(2,1) - m(3) * pkin(1) - t19 * mrSges(3,1) + sin(pkin(8)) * mrSges(3,2);
t17 = pkin(8) + qJ(3);
t15 = sin(t17);
t16 = cos(t17);
t55 = t16 * pkin(3) + t15 * pkin(7);
t54 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t29 = pkin(4) * t23 + qJ(5) * t21;
t53 = (mrSges(4,2) + t60) * t16 + (-m(6) * (-pkin(3) - t29) + m(5) * pkin(3) + mrSges(4,1) + t59) * t15;
t33 = t16 * mrSges(4,1) - t15 * mrSges(4,2);
t52 = t60 * t15 - t33;
t51 = pkin(7) * t56;
t50 = m(6) * pkin(4) + t62;
t49 = m(6) * qJ(5) + t61;
t48 = g(3) * t15;
t24 = cos(qJ(1));
t45 = t15 * t24;
t44 = t16 * t24;
t22 = sin(qJ(1));
t43 = t22 * t21;
t42 = t22 * t23;
t40 = t24 * t21;
t39 = t24 * t23;
t12 = t19 * pkin(2) + pkin(1);
t20 = -pkin(6) - qJ(2);
t36 = t24 * t12 - t22 * t20;
t4 = t16 * t39 + t43;
t3 = t16 * t40 - t42;
t2 = t16 * t42 - t40;
t1 = t16 * t43 + t39;
t5 = [(-m(4) * t36 + t60 * t45 - t56 * (pkin(3) * t44 + pkin(7) * t45 + t36) - t50 * t4 - t49 * t3 + (-t33 + t57) * t24 + t54 * t22) * g(2) + (t50 * t2 + t49 * t1 + (m(4) * t12 - t56 * (-t12 - t55) - t52 - t57) * t22 + (t58 * t20 + t54) * t24) * g(1), (-g(1) * t22 + g(2) * t24) * (m(3) + t58), (-t56 * t55 + (-m(6) * t29 - t59) * t16 + t52) * g(3) + (t53 * t24 - t44 * t51) * g(1) + (-t16 * t51 + t53) * g(2) * t22, (t21 * t50 - t49 * t23) * t48 + (t1 * t50 - t49 * t2) * g(2) + (t3 * t50 - t49 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - t21 * t48) * m(6)];
taug = t5(:);
