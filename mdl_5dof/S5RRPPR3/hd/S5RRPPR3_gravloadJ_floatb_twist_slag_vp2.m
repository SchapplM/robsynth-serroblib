% Calculate Gravitation load on the joints for
% S5RRPPR3
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:29
% DurationCPUTime: 0.24s
% Computational Cost: add. (220->49), mult. (144->50), div. (0->0), fcn. (102->8), ass. (0->29)
t54 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t52 = mrSges(6,1) * t23 + mrSges(6,2) * t25;
t53 = -t52 + mrSges(4,2) - mrSges(5,3);
t51 = m(5) + m(6);
t22 = qJ(1) + qJ(2);
t18 = pkin(8) + t22;
t15 = sin(t18);
t16 = cos(t18);
t50 = -g(1) * t15 + g(2) * t16;
t19 = sin(t22);
t20 = cos(t22);
t49 = mrSges(3,1) * t19 + mrSges(3,2) * t20 + t53 * t16 + (-m(6) * (-pkin(3) - pkin(7)) - t54) * t15;
t48 = -t20 * mrSges(3,1) + t19 * mrSges(3,2) + t53 * t15 + t54 * t16;
t24 = sin(qJ(1));
t47 = pkin(1) * t24;
t46 = pkin(2) * t19;
t17 = pkin(2) * t20;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t39 = t16 * pkin(3) + t15 * qJ(4) + t17;
t6 = t16 * qJ(4);
t38 = t6 - t46;
t34 = t21 + t39;
t33 = -t46 - t47;
t28 = -pkin(3) * t15 + t38;
t12 = t16 * pkin(7);
t1 = [(-mrSges(2,1) * t26 + t24 * mrSges(2,2) - m(3) * t21 - m(4) * (t17 + t21) - m(5) * t34 - m(6) * (t12 + t34) + t48) * g(2) + (t24 * mrSges(2,1) + mrSges(2,2) * t26 + m(3) * t47 - m(4) * t33 - m(5) * (t28 - t47) - m(6) * (t33 + t6) + t49) * g(1), (-m(4) * t17 - m(5) * t39 - m(6) * (t12 + t39) + t48) * g(2) + (m(4) * t46 - m(5) * t28 - m(6) * t38 + t49) * g(1), (-m(4) - t51) * g(3), t51 * t50, g(3) * t52 + t50 * (mrSges(6,1) * t25 - mrSges(6,2) * t23)];
taug = t1(:);
