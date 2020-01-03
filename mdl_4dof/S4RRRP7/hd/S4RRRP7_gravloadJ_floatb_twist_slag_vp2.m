% Calculate Gravitation load on the joints for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:09
% EndTime: 2019-12-31 17:20:11
% DurationCPUTime: 0.42s
% Computational Cost: add. (131->56), mult. (297->71), div. (0->0), fcn. (282->6), ass. (0->30)
t50 = mrSges(4,1) + mrSges(5,1);
t49 = -mrSges(4,2) + mrSges(5,3);
t48 = -mrSges(4,3) - mrSges(5,2);
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t47 = t49 * t16 + t50 * t19;
t46 = -m(4) - m(5);
t45 = mrSges(2,2) - mrSges(3,3);
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t44 = t20 * pkin(2) + t17 * pkin(6);
t25 = pkin(3) * t19 + qJ(4) * t16;
t43 = (mrSges(3,2) + t48) * t20 + (-m(5) * (-pkin(2) - t25) + m(4) * pkin(2) + mrSges(3,1) + t47) * t17;
t29 = t20 * mrSges(3,1) - t17 * mrSges(3,2);
t42 = t48 * t17 - t29;
t41 = pkin(6) * t46;
t40 = m(5) * pkin(3) + t50;
t39 = m(5) * qJ(4) + t49;
t38 = g(3) * t17;
t18 = sin(qJ(1));
t35 = t18 * t20;
t21 = cos(qJ(1));
t34 = t21 * t17;
t33 = t21 * t20;
t31 = t21 * pkin(1) + t18 * pkin(5);
t4 = t18 * t16 + t19 * t33;
t3 = t16 * t33 - t18 * t19;
t2 = -t21 * t16 + t19 * t35;
t1 = t16 * t35 + t21 * t19;
t5 = [(-m(3) * t31 + t46 * (pkin(2) * t33 + pkin(6) * t34 + t31) - t40 * t4 + t48 * t34 - t39 * t3 + (-mrSges(2,1) - t29) * t21 + t45 * t18) * g(2) + (t40 * t2 + t39 * t1 + (m(3) * pkin(1) + mrSges(2,1) + t46 * (-pkin(1) - t44) - t42) * t18 + (t45 + (-m(3) + t46) * pkin(5)) * t21) * g(1), (t46 * t44 + (-m(5) * t25 - t47) * t20 + t42) * g(3) + (t43 * t18 + t35 * t41) * g(2) + (t43 * t21 + t33 * t41) * g(1), (t40 * t16 - t39 * t19) * t38 + (t40 * t1 - t39 * t2) * g(2) + (t40 * t3 - t39 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - t16 * t38) * m(5)];
taug = t5(:);
