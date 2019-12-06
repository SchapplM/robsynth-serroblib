% Calculate Gravitation load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:06
% DurationCPUTime: 0.39s
% Computational Cost: add. (161->43), mult. (233->54), div. (0->0), fcn. (203->8), ass. (0->27)
t51 = m(5) + m(6);
t29 = m(4) + t51;
t55 = mrSges(5,1) + mrSges(6,1);
t54 = -mrSges(5,2) + mrSges(6,3);
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t53 = t54 * t17 + t55 * t19 + mrSges(4,1);
t52 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t40 = m(6) * pkin(4) + t55;
t39 = -m(6) * qJ(5) - t54;
t14 = qJ(2) + pkin(8);
t11 = sin(t14);
t34 = g(3) * t11;
t20 = cos(qJ(2));
t13 = t20 * pkin(2);
t15 = sin(pkin(7));
t33 = t15 * t17;
t32 = t15 * t19;
t16 = cos(pkin(7));
t31 = t16 * t17;
t30 = t16 * t19;
t24 = pkin(4) * t19 + qJ(5) * t17;
t18 = sin(qJ(2));
t12 = cos(t14);
t3 = t12 * t31 - t32;
t1 = t12 * t33 + t30;
t2 = [(-m(2) - m(3) - t29) * g(3), (-m(4) * t13 - mrSges(3,1) * t20 + t18 * mrSges(3,2) - t51 * (t12 * pkin(3) + t11 * pkin(6) + t13) + (-m(6) * t24 - t53) * t12 + t52 * t11) * g(3) + (mrSges(3,2) * t20 + (-m(6) * (-pkin(3) - t24) + m(5) * pkin(3) + t53) * t11 + (-pkin(6) * t51 + t52) * t12 + (t29 * pkin(2) + mrSges(3,1)) * t18) * (g(1) * t16 + g(2) * t15), (-g(1) * t15 + g(2) * t16) * t29, (t40 * t17 + t39 * t19) * t34 + (t39 * (t12 * t32 - t31) + t40 * t1) * g(2) + (t39 * (t12 * t30 + t33) + t40 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - t17 * t34) * m(6)];
taug = t2(:);
