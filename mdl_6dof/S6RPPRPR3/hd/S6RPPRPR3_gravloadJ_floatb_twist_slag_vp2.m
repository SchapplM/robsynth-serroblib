% Calculate Gravitation load on the joints for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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

function taug = S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:40:19
% EndTime: 2018-11-23 15:40:20
% DurationCPUTime: 0.39s
% Computational Cost: add. (275->69), mult. (250->80), div. (0->0), fcn. (201->10), ass. (0->38)
t59 = -mrSges(6,2) + mrSges(7,3);
t55 = m(6) + m(7);
t16 = qJ(4) + pkin(10);
t11 = sin(t16);
t13 = cos(t16);
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t58 = -t20 * mrSges(5,1) - t11 * mrSges(6,1) - t23 * mrSges(5,2) + t59 * t13;
t17 = qJ(1) + pkin(9);
t12 = sin(t17);
t14 = cos(t17);
t57 = -g(1) * t12 + g(2) * t14;
t56 = -m(4) - m(5);
t53 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t47 = t13 * pkin(8);
t52 = mrSges(3,2) - mrSges(4,3) - m(7) * (t11 * pkin(5) - t47) + t58;
t46 = t20 * pkin(4);
t21 = sin(qJ(1));
t45 = t21 * pkin(1);
t24 = cos(qJ(1));
t15 = t24 * pkin(1);
t19 = sin(qJ(6));
t43 = t12 * t19;
t22 = cos(qJ(6));
t42 = t12 * t22;
t40 = t14 * t19;
t39 = t14 * t22;
t38 = t14 * pkin(2) + t12 * qJ(3) + t15;
t37 = t55 - t56;
t36 = t14 * qJ(3) - t45;
t30 = -t12 * pkin(2) + t36;
t27 = m(7) * pkin(5) + t22 * mrSges(7,1) - t19 * mrSges(7,2);
t18 = -qJ(5) - pkin(7);
t4 = t11 * t39 - t43;
t3 = t11 * t40 + t42;
t2 = t11 * t42 + t40;
t1 = -t11 * t43 + t39;
t5 = [(-m(3) * t15 - t24 * mrSges(2,1) - t2 * mrSges(7,1) + t21 * mrSges(2,2) - t1 * mrSges(7,2) + t56 * t38 - t55 * (t12 * t46 - t14 * t18 + t38) + (-m(5) * pkin(7) - t53) * t14 + t52 * t12) * g(2) + (m(3) * t45 - m(4) * t30 - m(5) * t36 + t21 * mrSges(2,1) - t4 * mrSges(7,1) + t24 * mrSges(2,2) + t3 * mrSges(7,2) - t55 * (t12 * t18 + t14 * t46 + t30) + (-m(5) * (-pkin(2) - pkin(7)) + t53) * t12 + t52 * t14) * g(1) (-m(3) - t37) * g(3), t57 * t37 (m(6) * t46 - m(7) * (-t46 + t47) + t27 * t11 - t58) * g(3) + t57 * (-mrSges(5,2) * t20 + (mrSges(6,1) + t27) * t13 + (m(7) * pkin(8) + t59) * t11 + (t55 * pkin(4) + mrSges(5,1)) * t23) t55 * (-g(1) * t14 - g(2) * t12) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - g(3) * (-mrSges(7,1) * t19 - mrSges(7,2) * t22) * t13];
taug  = t5(:);
