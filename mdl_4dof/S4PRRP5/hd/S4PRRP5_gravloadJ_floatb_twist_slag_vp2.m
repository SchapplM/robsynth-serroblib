% Calculate Gravitation load on the joints for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:52
% DurationCPUTime: 0.23s
% Computational Cost: add. (70->28), mult. (163->38), div. (0->0), fcn. (140->6), ass. (0->15)
t32 = m(5) + m(4);
t28 = mrSges(4,2) + mrSges(5,2);
t27 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t31 = t32 * pkin(2) + t27 * t11 - t28 * t9 + mrSges(3,1);
t30 = mrSges(3,2) + m(5) * (-qJ(4) - pkin(5)) - mrSges(5,3) - m(4) * pkin(5) - mrSges(4,3);
t6 = sin(pkin(6));
t7 = cos(pkin(6));
t29 = g(1) * t7 + g(2) * t6;
t12 = cos(qJ(2));
t21 = t12 * t9;
t20 = t12 * t11;
t10 = sin(qJ(2));
t1 = [(-m(2) - m(3) - t32) * g(3), (-t31 * g(3) + t29 * t30) * t12 + (t30 * g(3) + t29 * t31) * t10, (t28 * t11 + t27 * t9) * g(3) * t10 + (-t28 * (-t6 * t20 + t7 * t9) - t27 * (-t7 * t11 - t6 * t21)) * g(2) + (-t28 * (-t7 * t20 - t6 * t9) - t27 * (t6 * t11 - t7 * t21)) * g(1), (g(3) * t12 - t29 * t10) * m(5)];
taug = t1(:);
