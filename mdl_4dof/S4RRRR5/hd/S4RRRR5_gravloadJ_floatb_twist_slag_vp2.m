% Calculate Gravitation load on the joints for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:22
% EndTime: 2019-12-31 17:27:24
% DurationCPUTime: 0.40s
% Computational Cost: add. (163->70), mult. (276->86), div. (0->0), fcn. (255->8), ass. (0->41)
t60 = mrSges(4,3) + mrSges(5,3);
t23 = cos(qJ(3));
t13 = t23 * pkin(3) + pkin(2);
t19 = qJ(3) + qJ(4);
t14 = sin(t19);
t15 = cos(t19);
t20 = sin(qJ(3));
t59 = -m(4) * pkin(2) - m(5) * t13 - t23 * mrSges(4,1) - t15 * mrSges(5,1) + t20 * mrSges(4,2) + t14 * mrSges(5,2);
t26 = -pkin(7) - pkin(6);
t58 = -m(4) * pkin(6) + m(5) * t26 - t60;
t52 = m(5) * pkin(3);
t57 = mrSges(4,1) + t52;
t56 = mrSges(2,2) - mrSges(3,3);
t55 = -m(3) - m(4) - m(5);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t34 = t24 * mrSges(3,1) - t21 * mrSges(3,2);
t53 = t60 * t21 + mrSges(2,1) + t34;
t25 = cos(qJ(1));
t41 = t25 * t15;
t22 = sin(qJ(1));
t43 = t22 * t24;
t5 = t14 * t43 + t41;
t42 = t25 * t14;
t6 = -t15 * t43 + t42;
t51 = -t5 * mrSges(5,1) + t6 * mrSges(5,2);
t7 = t22 * t15 - t24 * t42;
t8 = t22 * t14 + t24 * t41;
t50 = t7 * mrSges(5,1) - t8 * mrSges(5,2);
t47 = g(3) * t21;
t44 = t22 * t20;
t40 = t25 * t20;
t39 = t25 * t23;
t35 = t24 * pkin(2) + t21 * pkin(6);
t32 = -mrSges(5,1) * t14 - mrSges(5,2) * t15;
t31 = t24 * t13 - t21 * t26;
t11 = t22 * t23 - t24 * t40;
t9 = t20 * t43 + t39;
t12 = t24 * t39 + t44;
t10 = -t23 * t43 + t40;
t1 = [(-t44 * t52 - t12 * mrSges(4,1) - t8 * mrSges(5,1) - t11 * mrSges(4,2) - t7 * mrSges(5,2) + t55 * (t25 * pkin(1) + t22 * pkin(5)) + t56 * t22 + (-m(4) * t35 - m(5) * t31 - t53) * t25) * g(2) + (-t40 * t52 - t10 * mrSges(4,1) - t6 * mrSges(5,1) - t9 * mrSges(4,2) - t5 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t35) - m(5) * (-pkin(1) - t31) + t53) * t22 + (t55 * pkin(5) + t56) * t25) * g(1), (t58 * t21 + t59 * t24 - t34) * g(3) + (g(1) * t25 + g(2) * t22) * ((mrSges(3,2) + t58) * t24 + (mrSges(3,1) - t59) * t21), (mrSges(4,2) * t23 + t57 * t20 - t32) * t47 + (-mrSges(4,2) * t10 + t57 * t9 - t51) * g(2) + (mrSges(4,2) * t12 - t57 * t11 - t50) * g(1), -g(1) * t50 - g(2) * t51 - t32 * t47];
taug = t1(:);
