% Calculate Gravitation load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:04:48
% DurationCPUTime: 0.35s
% Computational Cost: add. (130->43), mult. (211->71), div. (0->0), fcn. (202->10), ass. (0->27)
t35 = m(5) + m(6);
t40 = pkin(3) * t35 + mrSges(4,1);
t14 = sin(qJ(5));
t16 = cos(qJ(5));
t37 = m(6) * pkin(4) + mrSges(6,1) * t16 - mrSges(6,2) * t14 + mrSges(5,1);
t36 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t10 = sin(pkin(8));
t32 = t10 * t14;
t31 = t10 * t16;
t11 = sin(pkin(7));
t12 = cos(pkin(8));
t30 = t11 * t12;
t15 = sin(qJ(3));
t29 = t11 * t15;
t17 = cos(qJ(3));
t28 = t11 * t17;
t13 = cos(pkin(7));
t27 = t12 * t13;
t26 = t13 * t15;
t25 = t13 * t17;
t22 = m(3) + m(4) + t35;
t9 = qJ(3) + pkin(9);
t8 = cos(t9);
t7 = sin(t9);
t4 = t11 * t7 + t27 * t8;
t2 = -t13 * t7 + t30 * t8;
t1 = [(-m(2) - t22) * g(3), (-g(1) * t11 + g(2) * t13) * t22, (mrSges(4,2) * t17 + t15 * t40 + t36 * t8 + t37 * t7) * g(3) * t10 + (-(-t12 * t25 - t29) * mrSges(4,2) + t36 * t4 - t37 * (t11 * t8 - t27 * t7) + t40 * (t12 * t26 - t28)) * g(1) + (-(-t12 * t28 + t26) * mrSges(4,2) - t37 * (-t13 * t8 - t30 * t7) + t36 * t2 - t40 * (-t12 * t29 - t25)) * g(2), (g(3) * t12 + t10 * (-g(1) * t13 - g(2) * t11)) * t35, -g(1) * ((t13 * t31 - t14 * t4) * mrSges(6,1) + (-t13 * t32 - t16 * t4) * mrSges(6,2)) - g(2) * ((t11 * t31 - t14 * t2) * mrSges(6,1) + (-t11 * t32 - t16 * t2) * mrSges(6,2)) - g(3) * ((-t12 * t16 - t32 * t8) * mrSges(6,1) + (t12 * t14 - t31 * t8) * mrSges(6,2))];
taug = t1(:);
