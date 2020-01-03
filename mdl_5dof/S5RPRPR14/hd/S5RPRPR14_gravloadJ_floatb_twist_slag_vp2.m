% Calculate Gravitation load on the joints for
% S5RPRPR14
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (155->59), mult. (227->71), div. (0->0), fcn. (187->8), ass. (0->32)
t53 = -mrSges(5,2) + mrSges(6,3);
t49 = m(5) + m(6);
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t13 = qJ(3) + pkin(8);
t8 = sin(t13);
t9 = cos(t13);
t52 = -t16 * mrSges(4,1) - t8 * mrSges(5,1) - t19 * mrSges(4,2) + t53 * t9;
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t51 = -g(1) * t17 + g(2) * t20;
t50 = -m(3) - m(4);
t47 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t44 = pkin(7) * t9;
t46 = mrSges(2,2) - mrSges(3,3) - m(6) * (pkin(4) * t8 - t44) + t52;
t43 = pkin(3) * t16;
t15 = sin(qJ(5));
t37 = t15 * t20;
t36 = t17 * t15;
t18 = cos(qJ(5));
t35 = t17 * t18;
t34 = t18 * t20;
t33 = t20 * pkin(1) + t17 * qJ(2);
t11 = t20 * qJ(2);
t31 = -t17 * pkin(1) + t11;
t23 = m(6) * pkin(4) + t18 * mrSges(6,1) - t15 * mrSges(6,2);
t14 = -qJ(4) - pkin(6);
t4 = t8 * t34 - t36;
t3 = t37 * t8 + t35;
t2 = t8 * t35 + t37;
t1 = -t8 * t36 + t34;
t5 = [(-t2 * mrSges(6,1) - t1 * mrSges(6,2) + t50 * t33 - t49 * (-t14 * t20 + t17 * t43 + t33) + (-m(4) * pkin(6) - t47) * t20 + t46 * t17) * g(2) + (-m(3) * t31 - m(4) * t11 - t4 * mrSges(6,1) + t3 * mrSges(6,2) - t49 * (t17 * t14 + t20 * t43 + t31) + (-m(4) * (-pkin(1) - pkin(6)) + t47) * t17 + t46 * t20) * g(1), t51 * (t49 - t50), (m(5) * t43 - m(6) * (-t43 + t44) + t23 * t8 - t52) * g(3) + t51 * (-mrSges(4,2) * t16 + (mrSges(5,1) + t23) * t9 + (m(6) * pkin(7) + t53) * t8 + (t49 * pkin(3) + mrSges(4,1)) * t19), t49 * (-g(1) * t20 - g(2) * t17), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(3) * (-mrSges(6,1) * t15 - mrSges(6,2) * t18) * t9];
taug = t5(:);
