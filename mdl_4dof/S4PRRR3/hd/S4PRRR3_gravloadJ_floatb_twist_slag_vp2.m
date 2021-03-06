% Calculate Gravitation load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->22), mult. (77->24), div. (0->0), fcn. (54->6), ass. (0->16)
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t26 = t16 * mrSges(5,1) - t15 * mrSges(5,2);
t25 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t24 = m(5) * pkin(3) + mrSges(4,1) + t26;
t22 = m(4) + m(5);
t23 = t22 * pkin(2) + mrSges(3,1);
t14 = pkin(7) + qJ(2);
t13 = qJ(3) + t14;
t10 = cos(t13);
t9 = sin(t13);
t18 = -t24 * t10 + t25 * t9;
t17 = t25 * t10 + t24 * t9;
t12 = cos(t14);
t11 = sin(t14);
t1 = [(-m(2) - m(3) - t22) * g(3), (mrSges(3,2) * t11 - t23 * t12 + t18) * g(2) + (mrSges(3,2) * t12 + t23 * t11 + t17) * g(1), t17 * g(1) + t18 * g(2), -g(3) * t26 + (g(1) * t10 + g(2) * t9) * (mrSges(5,1) * t15 + mrSges(5,2) * t16)];
taug = t1(:);
