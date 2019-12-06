% Calculate Gravitation load on the joints for
% S5PRPRP1
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:18
% EndTime: 2019-12-05 15:28:19
% DurationCPUTime: 0.19s
% Computational Cost: add. (164->41), mult. (136->39), div. (0->0), fcn. (99->6), ass. (0->19)
t27 = m(5) + m(6);
t8 = pkin(7) + qJ(2);
t4 = sin(t8);
t6 = cos(t8);
t26 = g(1) * t6 + g(2) * t4;
t10 = cos(pkin(8));
t7 = pkin(8) + qJ(4);
t3 = sin(t7);
t5 = cos(t7);
t18 = mrSges(5,1) * t5 - mrSges(5,2) * t3;
t25 = -mrSges(3,1) - m(4) * pkin(2) - t10 * mrSges(4,1) + sin(pkin(8)) * mrSges(4,2) - t18;
t24 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t21 = m(4) + t27;
t16 = t5 * mrSges(6,1) + t3 * mrSges(6,3);
t15 = pkin(4) * t5 + qJ(5) * t3;
t13 = m(6) * t15 + t16;
t11 = -pkin(6) - qJ(3);
t2 = pkin(3) * t10 + pkin(2);
t1 = [(-m(2) - m(3) - t21) * g(3), (-t27 * (-t11 * t4 + t6 * t2) + (-t13 + t25) * t6 + t24 * t4) * g(2) + ((t27 * t11 + t24) * t6 + (m(5) * t2 - m(6) * (-t15 - t2) + t16 - t25) * t4) * g(1), (-g(1) * t4 + g(2) * t6) * t21, (-t13 - t18) * g(3) + ((-m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3)) * t5 + (m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1)) * t3) * t26, (g(3) * t5 - t26 * t3) * m(6)];
taug = t1(:);
