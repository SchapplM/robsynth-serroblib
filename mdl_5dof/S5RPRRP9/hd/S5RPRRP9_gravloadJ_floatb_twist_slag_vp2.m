% Calculate Gravitation load on the joints for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:36
% EndTime: 2019-12-31 18:48:37
% DurationCPUTime: 0.34s
% Computational Cost: add. (244->61), mult. (224->63), div. (0->0), fcn. (171->8), ass. (0->31)
t62 = -mrSges(5,1) - mrSges(6,1);
t21 = pkin(8) + qJ(3);
t16 = sin(t21);
t17 = cos(t21);
t18 = qJ(4) + t21;
t13 = sin(t18);
t14 = cos(t18);
t55 = t62 * t14 + (mrSges(5,2) - mrSges(6,3)) * t13;
t61 = -mrSges(4,1) * t17 + mrSges(4,2) * t16 + t55;
t56 = t14 * pkin(4) + t13 * qJ(5);
t60 = m(6) * t56;
t59 = (-m(6) * qJ(5) - mrSges(6,3)) * t14;
t58 = m(5) + m(6);
t25 = sin(qJ(1));
t26 = cos(qJ(1));
t54 = g(1) * t26 + g(2) * t25;
t23 = cos(pkin(8));
t15 = t23 * pkin(2) + pkin(1);
t53 = m(4) * t15 + mrSges(2,1) + m(3) * pkin(1) + t23 * mrSges(3,1) - sin(pkin(8)) * mrSges(3,2) - t61;
t24 = -pkin(6) - qJ(2);
t52 = -m(3) * qJ(2) + m(4) * t24 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t51 = pkin(3) * t16;
t12 = pkin(3) * t17;
t46 = mrSges(5,2) * t14;
t43 = t59 * t25;
t42 = t59 * t26;
t28 = t46 + (m(6) * pkin(4) - t62) * t13;
t27 = m(6) * (-pkin(4) * t13 - t51) - t13 * mrSges(6,1);
t20 = -pkin(7) + t24;
t2 = t12 + t15;
t1 = [(-t58 * (t26 * t2 - t25 * t20) + (-t53 - t60) * t26 + t52 * t25) * g(2) + ((t58 * t20 + t52) * t26 + (m(5) * t2 - m(6) * (-t2 - t56) + t53) * t25) * g(1), (-g(1) * t25 + g(2) * t26) * (m(3) + m(4) + t58), -g(1) * (t26 * t27 - t42) - g(2) * (t25 * t27 - t43) + (-m(5) * t12 - m(6) * (t12 + t56) + t61) * g(3) + (m(5) * t51 + mrSges(4,1) * t16 + mrSges(5,1) * t13 + mrSges(4,2) * t17 + t46) * t54, (t55 - t60) * g(3) + (t25 * t28 + t43) * g(2) + (t26 * t28 + t42) * g(1), (g(3) * t14 - t54 * t13) * m(6)];
taug = t1(:);
