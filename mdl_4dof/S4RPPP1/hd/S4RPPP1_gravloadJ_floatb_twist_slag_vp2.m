% Calculate Gravitation load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RPPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:25
% DurationCPUTime: 0.17s
% Computational Cost: add. (181->43), mult. (234->52), div. (0->0), fcn. (199->10), ass. (0->25)
t37 = m(4) + m(5);
t38 = t37 * qJ(3) - mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t23 = sin(qJ(1));
t24 = cos(qJ(1));
t20 = sin(pkin(4));
t35 = qJ(2) * t20;
t36 = t24 * pkin(1) + t23 * t35;
t18 = pkin(4) - pkin(6);
t13 = cos(t18) / 0.2e1;
t17 = pkin(4) + pkin(6);
t15 = cos(t17);
t33 = t13 + t15 / 0.2e1;
t12 = -sin(t18) / 0.2e1;
t14 = sin(t17);
t32 = t14 / 0.2e1 + t12;
t29 = -t23 * pkin(1) + t24 * t35;
t27 = m(5) * qJ(4) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t25 = mrSges(2,2) + (-m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1) - mrSges(3,3)) * t20;
t21 = cos(pkin(6));
t19 = sin(pkin(6));
t6 = t21 * t24 - t23 * t32;
t5 = t24 * t19 + t23 * t33;
t4 = t23 * t21 + t24 * t32;
t3 = t19 * t23 - t24 * t33;
t1 = [(-mrSges(2,1) * t24 - m(3) * t36 - t27 * t6 - t38 * t5 + t25 * t23 - t37 * (t6 * pkin(2) + t36)) * g(2) + (t23 * mrSges(2,1) - m(3) * t29 + t27 * t4 + t38 * t3 + t25 * t24 + t37 * (t4 * pkin(2) - t29)) * g(1) (-cos(pkin(4)) * g(3) + (-g(1) * t23 + g(2) * t24) * t20) * (m(3) + t37) t37 * (-g(1) * t5 - g(2) * t3 - g(3) * (-t14 / 0.2e1 + t12)) (-g(1) * t6 - g(2) * t4 - g(3) * (t13 - t15 / 0.2e1)) * m(5)];
taug  = t1(:);
