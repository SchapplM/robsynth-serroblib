% Calculate Gravitation load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:54
% EndTime: 2019-12-31 18:25:55
% DurationCPUTime: 0.32s
% Computational Cost: add. (197->60), mult. (241->64), div. (0->0), fcn. (240->8), ass. (0->38)
t19 = sin(qJ(5));
t22 = cos(qJ(5));
t27 = -mrSges(6,1) * t22 + mrSges(6,2) * t19;
t55 = m(6) * pkin(4) + mrSges(5,1) - t27;
t54 = -mrSges(5,2) + mrSges(6,3);
t53 = -m(3) - m(4);
t52 = mrSges(2,1) + mrSges(3,1);
t51 = mrSges(2,2) - mrSges(3,3);
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t39 = qJ(3) + pkin(8);
t33 = sin(t39);
t34 = cos(t39);
t1 = -t21 * t33 - t24 * t34;
t2 = -t21 * t34 + t24 * t33;
t23 = cos(qJ(3));
t41 = t24 * t23;
t20 = sin(qJ(3));
t43 = t21 * t20;
t3 = -t41 - t43;
t42 = t21 * t23;
t44 = t20 * t24;
t4 = -t42 + t44;
t50 = -t3 * mrSges(4,1) - t4 * mrSges(4,2) - t55 * t1 + t54 * t2;
t49 = t4 * mrSges(4,1) - t3 * mrSges(4,2) + t54 * t1 + t55 * t2;
t48 = t2 * pkin(7);
t47 = m(5) + m(6);
t40 = t24 * pkin(1) + t21 * qJ(2);
t11 = pkin(3) * t44;
t37 = -pkin(7) * t1 - t11;
t17 = t24 * qJ(2);
t36 = -pkin(1) * t21 + t17;
t15 = pkin(3) * t23 + pkin(2);
t9 = pkin(3) * t43;
t35 = t24 * t15 + t40 + t9;
t32 = -pkin(3) * t41 - t9;
t10 = pkin(3) * t42;
t5 = [(-m(5) * t35 - m(6) * (t35 + t48) + t53 * t40 + (-m(4) * pkin(2) - t52) * t24 + t51 * t21 - t50) * g(2) + (-m(3) * t36 - m(4) * t17 - m(5) * (t11 + t17) - m(6) * (t36 - t37) + t51 * t24 + (-m(4) * (-pkin(1) - pkin(2)) - m(5) * (-pkin(1) - t15) + m(6) * t15 + t52) * t21 - t49) * g(1), (-g(1) * t21 + g(2) * t24) * (t47 - t53), (-m(5) * t32 - m(6) * (t32 - t48) + t50) * g(2) + (-m(5) * (-t11 + t10) - m(6) * (t10 + t37) + t49) * g(1), t47 * g(3), -g(3) * t27 + (-g(1) * t1 - g(2) * t2) * (mrSges(6,1) * t19 + mrSges(6,2) * t22)];
taug = t5(:);
