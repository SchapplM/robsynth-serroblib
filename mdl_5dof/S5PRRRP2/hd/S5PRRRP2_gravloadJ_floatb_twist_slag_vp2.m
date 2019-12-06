% Calculate Gravitation load on the joints for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (237->46), mult. (166->47), div. (0->0), fcn. (125->6), ass. (0->24)
t54 = mrSges(5,1) + mrSges(6,1);
t53 = -mrSges(5,2) + mrSges(6,3);
t22 = sin(qJ(4));
t23 = cos(qJ(4));
t52 = t53 * t22 + t54 * t23;
t51 = -mrSges(4,1) - t52;
t50 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t28 = pkin(4) * t23 + t22 * qJ(5);
t21 = pkin(8) + qJ(2);
t20 = qJ(3) + t21;
t16 = sin(t20);
t17 = cos(t20);
t48 = g(1) * t17 + g(2) * t16;
t47 = t50 * t17 + (-m(6) * (-pkin(3) - t28) - t51) * t16;
t46 = t50 * t16 + t51 * t17;
t18 = sin(t21);
t45 = pkin(2) * t18;
t19 = cos(t21);
t15 = pkin(2) * t19;
t37 = t17 * pkin(3) + t16 * pkin(7);
t13 = t17 * pkin(7);
t34 = -t16 * pkin(3) + t13;
t33 = t28 * t17 + t37;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (-mrSges(3,1) * t19 + mrSges(3,2) * t18 - m(4) * t15 - m(5) * (t15 + t37) - m(6) * (t15 + t33) + t46) * g(2) + (mrSges(3,1) * t18 + mrSges(3,2) * t19 + m(4) * t45 - m(5) * (t34 - t45) - m(6) * (t13 - t45) + t47) * g(1), (-m(5) * t37 - m(6) * t33 + t46) * g(2) + (-m(5) * t34 - m(6) * t13 + t47) * g(1), (-m(6) * t28 - t52) * g(3) + ((-m(6) * qJ(5) - t53) * t23 + (m(6) * pkin(4) + t54) * t22) * t48, (g(3) * t23 - t48 * t22) * m(6)];
taug = t1(:);
