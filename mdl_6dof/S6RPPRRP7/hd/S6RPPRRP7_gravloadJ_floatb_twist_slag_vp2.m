% Calculate Gravitation load on the joints for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2018-11-23 15:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:53
% EndTime: 2018-11-23 15:46:53
% DurationCPUTime: 0.51s
% Computational Cost: add. (254->64), mult. (347->73), div. (0->0), fcn. (296->8), ass. (0->35)
t61 = -mrSges(5,2) + m(6) * pkin(8) - m(7) * (-qJ(6) - pkin(8)) + mrSges(6,3) + mrSges(7,3);
t14 = pkin(9) + qJ(4);
t10 = cos(t14);
t9 = sin(t14);
t60 = -t9 * mrSges(5,1) + t10 * t61;
t59 = mrSges(6,1) + mrSges(7,1);
t54 = mrSges(6,2) + mrSges(7,2);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t51 = -g(1) * t20 + g(2) * t22;
t21 = cos(qJ(5));
t58 = -m(6) * pkin(4) - m(7) * (t21 * pkin(5) + pkin(4));
t46 = m(7) * pkin(5);
t19 = sin(qJ(5));
t53 = -t19 * t54 + t21 * t59 - t58;
t50 = -m(5) - m(6) - m(7);
t49 = t46 + t59;
t48 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t15 = sin(pkin(9));
t47 = -t15 * mrSges(4,1) - cos(pkin(9)) * mrSges(4,2) + mrSges(2,2) - mrSges(3,3) + t58 * t9 + t60;
t45 = pkin(3) * t15;
t41 = t20 * t19;
t40 = t20 * t21;
t39 = t22 * t19;
t38 = t22 * t21;
t37 = pkin(1) * t22 + qJ(2) * t20;
t36 = -m(4) + t50;
t12 = t22 * qJ(2);
t33 = -t20 * pkin(1) + t12;
t3 = t39 * t9 + t40;
t1 = -t41 * t9 + t38;
t18 = -pkin(7) - qJ(3);
t4 = t38 * t9 - t41;
t2 = t40 * t9 + t39;
t5 = [(-t39 * t46 + (-m(3) - m(4)) * t37 + t50 * (-t22 * t18 + t20 * t45 + t37) - t59 * t2 - t54 * t1 + (-m(4) * qJ(3) - t48) * t22 + t47 * t20) * g(2) + (t41 * t46 - m(3) * t33 - m(4) * t12 - t59 * t4 + t54 * t3 + t50 * (t18 * t20 + t22 * t45 + t33) + (-m(4) * (-pkin(1) - qJ(3)) + t48) * t20 + t47 * t22) * g(1), t51 * (m(3) - t36) (g(1) * t22 + g(2) * t20) * t36 (t53 * t9 - t60) * g(3) + t51 * (t61 * t9 + (mrSges(5,1) + t53) * t10) (t19 * t49 + t21 * t54) * g(3) * t10 + (-t3 * t49 - t4 * t54) * g(2) + (-t1 * t49 + t2 * t54) * g(1) (-g(3) * t9 - t10 * t51) * m(7)];
taug  = t5(:);
