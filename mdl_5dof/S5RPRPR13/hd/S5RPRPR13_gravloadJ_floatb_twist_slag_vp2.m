% Calculate Gravitation load on the joints for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:57
% EndTime: 2019-12-31 18:31:59
% DurationCPUTime: 0.42s
% Computational Cost: add. (199->70), mult. (252->74), div. (0->0), fcn. (210->8), ass. (0->37)
t59 = -mrSges(4,1) + mrSges(5,2);
t58 = mrSges(4,2) - mrSges(5,3);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t51 = g(1) * t22 + g(2) * t20;
t47 = m(5) + m(6);
t15 = pkin(8) + qJ(3);
t13 = sin(t15);
t10 = t13 * qJ(4);
t14 = cos(t15);
t42 = t14 * t22;
t54 = pkin(3) * t42 + t22 * t10;
t52 = t58 * t13 + t59 * t14;
t17 = cos(pkin(8));
t49 = -mrSges(2,1) - m(3) * pkin(1) - t17 * mrSges(3,1) + sin(pkin(8)) * mrSges(3,2) + t52;
t18 = -pkin(6) - qJ(2);
t48 = -m(3) * qJ(2) - m(6) * (pkin(4) - t18) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t44 = g(3) * t14;
t11 = t14 * pkin(3);
t19 = sin(qJ(5));
t41 = t20 * t19;
t21 = cos(qJ(5));
t40 = t20 * t21;
t39 = t22 * t19;
t38 = t22 * t21;
t37 = t11 + t10;
t12 = t17 * pkin(2) + pkin(1);
t8 = t22 * t12;
t34 = -t20 * t18 + t8;
t32 = -t12 - t10;
t31 = m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3);
t27 = t19 * mrSges(6,1) + t21 * mrSges(6,2);
t4 = -t13 * t41 + t38;
t3 = t13 * t40 + t39;
t2 = t13 * t39 + t40;
t1 = t13 * t38 - t41;
t5 = [(-m(4) * t34 - m(5) * (t34 + t54) - m(6) * (pkin(7) * t42 + t54 + t8) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - mrSges(6,3) * t42 + t49 * t22 + t48 * t20) * g(2) + (-t4 * mrSges(6,1) + t3 * mrSges(6,2) + ((m(4) + m(5)) * t18 + t48) * t22 + (m(4) * t12 - m(5) * (t32 - t11) - m(6) * t32 - t31 * t14 - t49) * t20) * g(1), (-g(1) * t20 + g(2) * t22) * (m(3) + m(4) + t47), (-m(5) * t37 - m(6) * (pkin(7) * t14 + t37) - t14 * mrSges(6,3) - t27 * t13 + t52) * g(3) + ((m(5) * pkin(3) - t31 - t59) * t13 + (-qJ(4) * t47 - t27 + t58) * t14) * t51, (-t13 * t51 + t44) * t47, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - (-mrSges(6,1) * t21 + mrSges(6,2) * t19) * t44];
taug = t5(:);
