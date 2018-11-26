% Calculate Gravitation load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2018-11-23 15:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:59
% EndTime: 2018-11-23 15:40:00
% DurationCPUTime: 0.48s
% Computational Cost: add. (321->81), mult. (275->85), div. (0->0), fcn. (224->10), ass. (0->45)
t66 = -mrSges(5,1) + mrSges(6,2);
t65 = mrSges(5,2) - mrSges(6,3);
t19 = qJ(1) + pkin(9);
t14 = sin(t19);
t16 = cos(t19);
t56 = g(1) * t16 + g(2) * t14;
t52 = m(6) + m(7);
t62 = -m(5) - t52;
t60 = m(3) + m(4);
t18 = pkin(10) + qJ(4);
t13 = sin(t18);
t10 = t13 * qJ(5);
t15 = cos(t18);
t44 = t15 * t16;
t59 = pkin(4) * t44 + t16 * t10;
t57 = t65 * t13 + t66 * t15;
t21 = cos(pkin(10));
t54 = -mrSges(3,1) - m(4) * pkin(2) - t21 * mrSges(4,1) + sin(pkin(10)) * mrSges(4,2) + t57;
t53 = -m(4) * qJ(3) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t49 = g(3) * t15;
t11 = t15 * pkin(4);
t24 = sin(qJ(1));
t48 = t24 * pkin(1);
t26 = cos(qJ(1));
t17 = t26 * pkin(1);
t12 = t21 * pkin(3) + pkin(2);
t47 = t16 * t12 + t17;
t23 = sin(qJ(6));
t46 = t14 * t23;
t25 = cos(qJ(6));
t45 = t14 * t25;
t43 = t16 * t23;
t42 = t16 * t25;
t41 = t11 + t10;
t39 = m(4) - t62;
t38 = -t12 - t10;
t22 = -pkin(7) - qJ(3);
t37 = -t14 * t22 + t47;
t36 = m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3);
t31 = t23 * mrSges(7,1) + t25 * mrSges(7,2);
t4 = -t13 * t46 + t42;
t3 = t13 * t45 + t43;
t2 = t13 * t43 + t45;
t1 = t13 * t42 - t46;
t5 = [(-t26 * mrSges(2,1) + t24 * mrSges(2,2) - m(5) * t37 - m(6) * (t37 + t59) - m(7) * (pkin(8) * t44 + t47 + t59) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t44 - t60 * t17 + t54 * t16 + (-m(7) * (pkin(5) - t22) + t53) * t14) * g(2) + (t24 * mrSges(2,1) - t4 * mrSges(7,1) + t26 * mrSges(2,2) + t3 * mrSges(7,2) + t60 * t48 + t62 * (-t16 * t22 - t48) + (-m(7) * pkin(5) + t53) * t16 + (m(5) * t12 - m(6) * (t38 - t11) - m(7) * t38 - t15 * t36 - t54) * t14) * g(1) (-m(3) - t39) * g(3) (-g(1) * t14 + g(2) * t16) * t39 (-m(6) * t41 - m(7) * (t15 * pkin(8) + t41) - t15 * mrSges(7,3) - t31 * t13 + t57) * g(3) + ((m(6) * pkin(4) - t36 - t66) * t13 + (-qJ(5) * t52 - t31 + t65) * t15) * t56 (-t13 * t56 + t49) * t52, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - (-mrSges(7,1) * t25 + mrSges(7,2) * t23) * t49];
taug  = t5(:);
