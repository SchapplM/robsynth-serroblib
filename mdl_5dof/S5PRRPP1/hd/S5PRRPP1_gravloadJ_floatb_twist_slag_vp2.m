% Calculate Gravitation load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:11
% EndTime: 2019-12-05 16:06:13
% DurationCPUTime: 0.24s
% Computational Cost: add. (180->43), mult. (158->42), div. (0->0), fcn. (117->6), ass. (0->20)
t31 = m(5) + m(6);
t35 = -mrSges(5,1) - mrSges(6,1);
t34 = mrSges(5,2) - mrSges(6,3);
t8 = pkin(7) + qJ(2);
t3 = sin(t8);
t5 = cos(t8);
t33 = g(1) * t5 + g(2) * t3;
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t9 = qJ(3) + pkin(8);
t4 = sin(t9);
t6 = cos(t9);
t32 = -mrSges(4,1) * t12 + t11 * mrSges(4,2) + t34 * t4 + t35 * t6;
t28 = m(4) * pkin(2) + mrSges(3,1) - t32;
t27 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t7 = t12 * pkin(3);
t18 = pkin(4) * t6 + qJ(5) * t4;
t10 = -qJ(4) - pkin(6);
t2 = t7 + pkin(2);
t1 = [(-m(2) - m(3) - m(4) - t31) * g(3), (-t31 * (-t10 * t3 + t5 * t2) + (-m(6) * t18 - t28) * t5 + t27 * t3) * g(2) + ((t31 * t10 + t27) * t5 + (m(5) * t2 - m(6) * (-t18 - t2) + t28) * t3) * g(1), (-m(5) * t7 - m(6) * (t18 + t7) + t32) * g(3) + t33 * (mrSges(4,2) * t12 + (-m(6) * qJ(5) + t34) * t6 + (m(6) * pkin(4) - t35) * t4 + (t31 * pkin(3) + mrSges(4,1)) * t11), t31 * (-g(1) * t3 + g(2) * t5), (g(3) * t6 - t33 * t4) * m(6)];
taug = t1(:);
