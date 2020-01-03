% Calculate Gravitation load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:22
% DurationCPUTime: 0.40s
% Computational Cost: add. (134->51), mult. (301->64), div. (0->0), fcn. (287->6), ass. (0->26)
t43 = -m(5) - m(6);
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t42 = -mrSges(2,1) + (-mrSges(3,1) - mrSges(4,1)) * t17 + (mrSges(3,2) - mrSges(4,3)) * t16;
t26 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t41 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t40 = mrSges(2,2) - mrSges(4,2) + mrSges(6,2) - mrSges(3,3) + mrSges(5,3);
t18 = sin(qJ(4));
t37 = cos(qJ(4));
t30 = t16 * t37;
t38 = -t17 * t18 + t30;
t20 = cos(qJ(1));
t35 = t20 * t17;
t19 = sin(qJ(1));
t34 = t20 * pkin(1) + t19 * qJ(2);
t33 = qJ(3) * t16;
t32 = m(4) - t43;
t28 = -pkin(1) - t33;
t27 = pkin(2) * t35 + t20 * t33 + t34;
t5 = t16 * t18 + t17 * t37;
t14 = t20 * qJ(2);
t4 = t5 * t20;
t3 = t18 * t35 - t20 * t30;
t2 = t5 * t19;
t1 = t38 * t19;
t6 = [(-m(3) * t34 - m(4) * t27 + t43 * (pkin(3) * t35 - pkin(6) * t19 + t27) - t26 * t4 - t41 * t3 + t42 * t20 + t40 * t19) * g(2) + (t43 * (((-pkin(2) - pkin(3)) * t17 + t28) * t19 - pkin(6) * t20 + t14) + t26 * t2 + (-m(3) - m(4)) * t14 - t41 * t1 + (m(3) * pkin(1) - m(4) * (-pkin(2) * t17 + t28) - t42) * t19 + t40 * t20) * g(1), (-g(1) * t19 + g(2) * t20) * (m(3) + t32), (g(3) * t17 + (-g(1) * t20 - g(2) * t19) * t16) * t32, (t26 * t5 - t38 * t41) * g(3) + (-t26 * t1 - t2 * t41) * g(2) + (t26 * t3 - t4 * t41) * g(1), (-g(1) * t3 + g(2) * t1 - g(3) * t5) * m(6)];
taug = t6(:);
