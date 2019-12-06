% Calculate Gravitation load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:22
% DurationCPUTime: 0.34s
% Computational Cost: add. (104->36), mult. (250->47), div. (0->0), fcn. (217->6), ass. (0->20)
t34 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t33 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t45 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t44 = -t13 * t34 + t15 * t33 + mrSges(3,2) - mrSges(4,3);
t38 = -m(5) - m(6);
t11 = sin(pkin(7));
t12 = cos(pkin(7));
t37 = g(1) * t12 + g(2) * t11;
t24 = m(4) - t38;
t16 = cos(qJ(2));
t29 = g(3) * t16;
t14 = sin(qJ(2));
t28 = pkin(2) * t16 + qJ(3) * t14;
t27 = t14 * t13;
t26 = t14 * t15;
t3 = t11 * t26 + t12 * t13;
t1 = t11 * t13 - t12 * t26;
t2 = [(-m(2) - m(3) - t24) * g(3), (-m(4) * t28 + t38 * (pkin(6) * t16 + t28) + t45 * t16 + t44 * t14) * g(3) + ((m(4) * pkin(2) + t38 * (-pkin(2) - pkin(6)) - t45) * t14 + (-qJ(3) * t24 + t44) * t16) * t37, (-t14 * t37 + t29) * t24, (t13 * t33 + t15 * t34) * t29 + (t33 * (-t11 * t27 + t12 * t15) - t34 * t3) * g(2) + (-t33 * (t11 * t15 + t12 * t27) + t34 * t1) * g(1), (-g(1) * t1 + g(2) * t3 - t15 * t29) * m(6)];
taug = t2(:);
