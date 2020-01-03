% Calculate Gravitation load on the joints for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:33
% EndTime: 2019-12-31 18:44:35
% DurationCPUTime: 0.47s
% Computational Cost: add. (260->62), mult. (317->79), div. (0->0), fcn. (294->8), ass. (0->37)
t58 = mrSges(5,1) + mrSges(6,1);
t57 = -mrSges(5,2) + mrSges(6,3);
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t56 = t57 * t20 + t58 * t23;
t55 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t54 = -m(5) - m(6);
t53 = mrSges(3,2) - mrSges(4,3);
t21 = sin(qJ(3));
t16 = t21 * pkin(7);
t24 = cos(qJ(3));
t52 = t24 * pkin(3) + t16;
t51 = -m(4) + t54;
t29 = pkin(4) * t23 + qJ(5) * t20;
t50 = t55 * t24 + (-m(6) * (-pkin(3) - t29) + m(5) * pkin(3) + mrSges(4,1) + t56) * t21;
t49 = -t24 * mrSges(4,1) + t55 * t21;
t48 = pkin(7) * t54;
t47 = -mrSges(3,1) + t49;
t46 = m(6) * pkin(4) + t58;
t45 = m(6) * qJ(5) + t57;
t44 = g(3) * t21;
t22 = sin(qJ(1));
t43 = t22 * pkin(1);
t25 = cos(qJ(1));
t18 = t25 * pkin(1);
t19 = qJ(1) + pkin(8);
t15 = cos(t19);
t42 = t15 * t24;
t39 = t24 * t20;
t38 = t24 * t23;
t14 = sin(t19);
t36 = t15 * pkin(2) + t14 * pkin(6) + t18;
t4 = t14 * t20 + t15 * t38;
t3 = -t14 * t23 + t15 * t39;
t2 = t14 * t38 - t15 * t20;
t1 = t14 * t39 + t15 * t23;
t5 = [(-m(3) * t18 - m(4) * t36 - t25 * mrSges(2,1) + t22 * mrSges(2,2) + t54 * (pkin(3) * t42 + t15 * t16 + t36) - t46 * t4 - t45 * t3 + t53 * t14 + t47 * t15) * g(2) + (m(3) * t43 + t22 * mrSges(2,1) + t25 * mrSges(2,2) + t51 * (t15 * pkin(6) - t43) + t46 * t2 + t53 * t15 + t45 * t1 + (m(4) * pkin(2) + t54 * (-pkin(2) - t52) - t47) * t14) * g(1), (-m(3) + t51) * g(3), (t54 * t52 + (-m(6) * t29 - t56) * t24 + t49) * g(3) + (t50 * t15 + t42 * t48) * g(1) + (t24 * t48 + t50) * g(2) * t14, (t46 * t20 - t45 * t23) * t44 + (t46 * t1 - t45 * t2) * g(2) + (t46 * t3 - t45 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - t20 * t44) * m(6)];
taug = t5(:);
