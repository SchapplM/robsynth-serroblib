% Calculate Gravitation load on the joints for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:03
% EndTime: 2019-12-31 16:36:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (162->48), mult. (412->81), div. (0->0), fcn. (474->10), ass. (0->31)
t47 = m(4) + m(5);
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t50 = m(5) * pkin(3) + t23 * mrSges(5,1) - t20 * mrSges(5,2) + mrSges(4,1);
t44 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3);
t49 = -t20 * mrSges(5,1) - t23 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t43 = -t44 * t21 + t50 * t24 + mrSges(3,1);
t48 = pkin(2) * t47 + t43;
t41 = -t47 * pkin(6) + t49;
t19 = sin(pkin(4));
t22 = sin(qJ(2));
t40 = t19 * t22;
t39 = t19 * t24;
t25 = cos(qJ(2));
t38 = t19 * t25;
t36 = cos(pkin(4));
t35 = cos(pkin(8));
t18 = sin(pkin(8));
t33 = t18 * t36;
t32 = t19 * t35;
t29 = t36 * t35;
t12 = t36 * t21 + t22 * t39;
t10 = -t22 * t33 + t35 * t25;
t9 = t35 * t22 + t25 * t33;
t8 = t18 * t25 + t22 * t29;
t7 = t18 * t22 - t25 * t29;
t4 = t18 * t19 * t21 + t10 * t24;
t2 = -t21 * t32 + t8 * t24;
t1 = [(-m(2) - m(3) - t47) * g(3), (-t47 * (pkin(2) * t38 + pkin(6) * t40) + (t49 * t22 - t43 * t25) * t19) * g(3) + (t41 * t8 + t48 * t7) * g(2) + (t41 * t10 + t48 * t9) * g(1), (t44 * t12 - t50 * (-t21 * t40 + t36 * t24)) * g(3) + (t44 * t2 - t50 * (-t8 * t21 - t24 * t32)) * g(2) + (t44 * t4 - t50 * (-t10 * t21 + t18 * t39)) * g(1), -g(1) * ((-t4 * t20 + t9 * t23) * mrSges(5,1) + (-t9 * t20 - t4 * t23) * mrSges(5,2)) - g(2) * ((-t2 * t20 + t7 * t23) * mrSges(5,1) + (-t2 * t23 - t7 * t20) * mrSges(5,2)) - g(3) * ((-t12 * t20 - t23 * t38) * mrSges(5,1) + (-t12 * t23 + t20 * t38) * mrSges(5,2))];
taug = t1(:);
