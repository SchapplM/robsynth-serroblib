% Calculate Gravitation load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (89->24), mult. (175->29), div. (0->0), fcn. (177->6), ass. (0->16)
t13 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t18 = -mrSges(5,2) - mrSges(6,2);
t8 = sin(qJ(4));
t9 = cos(qJ(4));
t27 = t13 * t9 + t18 * t8;
t26 = m(6) + m(5);
t23 = t26 * pkin(3) + mrSges(4,1) + t27;
t22 = -mrSges(4,2) - m(6) * (-qJ(5) - pkin(6)) + mrSges(6,3) + m(5) * pkin(6) + mrSges(5,3);
t21 = m(3) + m(4) + t26;
t20 = cos(qJ(3));
t19 = sin(qJ(3));
t17 = cos(pkin(7));
t16 = sin(pkin(7));
t2 = -t16 * t20 + t17 * t19;
t1 = -t16 * t19 - t17 * t20;
t3 = [(-m(2) - t21) * g(3), (-t16 * g(1) + t17 * g(2)) * t21, (-t23 * t1 + t22 * t2) * g(2) + (t22 * t1 + t23 * t2) * g(1), t27 * g(3) + (g(1) * t1 + g(2) * t2) * (-t13 * t8 + t18 * t9), (-g(1) * t2 + g(2) * t1) * m(6)];
taug = t3(:);
