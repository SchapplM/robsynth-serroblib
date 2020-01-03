% Calculate Gravitation load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:41
% DurationCPUTime: 0.17s
% Computational Cost: add. (37->20), mult. (89->34), div. (0->0), fcn. (71->6), ass. (0->13)
t3 = sin(qJ(4));
t5 = cos(qJ(4));
t18 = m(5) * pkin(3) + t5 * mrSges(5,1) - t3 * mrSges(5,2) + mrSges(4,1);
t17 = -m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3);
t1 = sin(pkin(6));
t2 = cos(pkin(6));
t16 = -g(1) * t1 + g(2) * t2;
t4 = sin(qJ(3));
t12 = t4 * t3;
t11 = t4 * t5;
t10 = m(3) + m(4) + m(5);
t6 = cos(qJ(3));
t7 = [(-m(2) - t10) * g(3), t16 * t10, (t17 * t6 + t18 * t4) * g(3) + t16 * (-t17 * t4 + t18 * t6), -g(1) * ((-t1 * t12 + t2 * t5) * mrSges(5,1) + (-t1 * t11 - t2 * t3) * mrSges(5,2)) - g(2) * ((t1 * t5 + t2 * t12) * mrSges(5,1) + (-t1 * t3 + t2 * t11) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t3 - mrSges(5,2) * t5) * t6];
taug = t7(:);
