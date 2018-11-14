% Calculate Gravitation load on the joints for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:26
% DurationCPUTime: 0.15s
% Computational Cost: add. (72->34), mult. (96->40), div. (0->0), fcn. (84->6), ass. (0->19)
t25 = m(4) + m(5);
t22 = m(3) + t25;
t23 = mrSges(2,1) + mrSges(3,1);
t21 = pkin(6) + qJ(4);
t20 = cos(t21);
t19 = sin(t21);
t14 = sin(qJ(1));
t15 = cos(qJ(1));
t1 = -t14 * t19 - t15 * t20;
t2 = -t14 * t20 + t15 * t19;
t18 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t17 = t1 * mrSges(5,1) + t2 * mrSges(5,2);
t12 = sin(pkin(6));
t16 = -m(5) * pkin(3) * t12 + mrSges(2,2) - mrSges(3,3);
t13 = cos(pkin(6));
t8 = pkin(3) * t13 + pkin(2);
t4 = t14 * t12 + t13 * t15;
t3 = t12 * t15 - t14 * t13;
t5 = [(-t3 * mrSges(4,1) - t4 * mrSges(4,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(2)) - m(5) * (-pkin(1) - t8) + t23) * t14 + t18 + (-t22 * qJ(2) + t16) * t15) * g(1) + (-t4 * mrSges(4,1) + t3 * mrSges(4,2) + (-m(4) * pkin(2) - m(5) * t8 - t23) * t15 + t16 * t14 + t17 - t22 * (t15 * pkin(1) + t14 * qJ(2))) * g(2) (-g(1) * t14 + g(2) * t15) * t22, t25 * g(3), -g(1) * t18 - g(2) * t17];
taug  = t5(:);
