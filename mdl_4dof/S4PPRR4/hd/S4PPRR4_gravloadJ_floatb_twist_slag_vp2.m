% Calculate Gravitation load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (64->22), mult. (89->34), div. (0->0), fcn. (71->6), ass. (0->15)
t6 = sin(qJ(4));
t7 = cos(qJ(4));
t21 = m(5) * pkin(3) + t7 * mrSges(5,1) - t6 * mrSges(5,2) + mrSges(4,1);
t20 = -m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3);
t4 = sin(pkin(6));
t16 = t4 * t6;
t15 = t4 * t7;
t5 = cos(pkin(6));
t14 = t5 * t6;
t13 = t5 * t7;
t12 = m(3) + m(4) + m(5);
t3 = pkin(7) + qJ(3);
t2 = cos(t3);
t1 = sin(t3);
t8 = [(-m(2) - t12) * g(3), (-g(1) * t4 + g(2) * t5) * t12, (t20 * t1 - t21 * t2) * g(3) + (t21 * t1 + t20 * t2) * (g(1) * t5 + g(2) * t4), -g(1) * ((-t2 * t14 + t15) * mrSges(5,1) + (-t2 * t13 - t16) * mrSges(5,2)) - g(2) * ((-t2 * t16 - t13) * mrSges(5,1) + (-t2 * t15 + t14) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t6 - mrSges(5,2) * t7) * t1];
taug = t8(:);
