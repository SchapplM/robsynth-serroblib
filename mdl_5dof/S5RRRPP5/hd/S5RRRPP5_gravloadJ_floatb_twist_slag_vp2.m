% Calculate Gravitation load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:35
% EndTime: 2019-12-31 20:57:37
% DurationCPUTime: 0.48s
% Computational Cost: add. (246->83), mult. (297->88), div. (0->0), fcn. (235->6), ass. (0->40)
t72 = -mrSges(6,2) - mrSges(5,3);
t22 = qJ(2) + qJ(3);
t19 = sin(t22);
t20 = cos(t22);
t71 = (-mrSges(4,1) - mrSges(5,1)) * t20 + (mrSges(4,2) + t72) * t19;
t10 = t19 * qJ(4);
t26 = cos(qJ(1));
t52 = t20 * t26;
t69 = pkin(3) * t52 + t26 * t10;
t48 = qJ(4) * t20;
t3 = t26 * t48;
t68 = -m(6) * t3 + t72 * t52;
t24 = sin(qJ(1));
t1 = t24 * t48;
t67 = t72 * t20 * t24 - m(6) * t1;
t40 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t32 = t40 * t19;
t23 = sin(qJ(2));
t59 = pkin(2) * t23;
t64 = m(6) * t59 - m(5) * (-pkin(3) * t19 - t59) + t19 * mrSges(5,1) - t32;
t63 = g(1) * t26 + g(2) * t24;
t62 = -t20 * mrSges(6,1) + t71;
t25 = cos(qJ(2));
t37 = t25 * mrSges(3,1) - t23 * mrSges(3,2);
t61 = -m(3) * pkin(1) - mrSges(2,1) - t37 + t71;
t27 = -pkin(7) - pkin(6);
t60 = mrSges(2,2) - m(3) * pkin(6) - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) - m(6) * (-qJ(5) - t27) + mrSges(6,3);
t17 = t20 * pkin(3);
t21 = t25 * pkin(2);
t55 = t19 * t24;
t54 = t19 * t26;
t49 = t17 + t10;
t47 = t21 + t49;
t18 = t21 + pkin(1);
t4 = t26 * t18;
t43 = -t24 * t27 + t4;
t41 = -t18 - t10;
t34 = mrSges(4,1) * t19 + mrSges(4,2) * t20;
t16 = t20 * pkin(4);
t2 = [(-m(4) * t43 - m(5) * (t43 + t69) - m(6) * (t4 + t69) + (-(m(6) * pkin(4) + mrSges(6,1)) * t20 + t61) * t26 + t60 * t24) * g(2) + (((m(4) + m(5)) * t27 + t60) * t26 + (m(4) * t18 - m(5) * (t41 - t17) - m(6) * t41 - t20 * t40 - t61) * t24) * g(1), (-m(5) * t1 + t64 * t24 + t67) * g(2) + (-m(5) * t3 + t64 * t26 + t68) * g(1) + (-t37 - m(4) * t21 - m(5) * t47 - m(6) * (t16 + t47) + t62) * g(3) + (m(4) * t59 + mrSges(3,1) * t23 + mrSges(3,2) * t25 + t34) * t63, t63 * t34 + (-m(5) * (-pkin(3) * t55 + t1) + mrSges(5,1) * t55 - t24 * t32 + t67) * g(2) + (-m(5) * (-pkin(3) * t54 + t3) + mrSges(5,1) * t54 - t26 * t32 + t68) * g(1) + (-m(5) * t49 - m(6) * (t16 + t49) + t62) * g(3), (g(3) * t20 - t19 * t63) * (m(5) + m(6)), (g(1) * t24 - g(2) * t26) * m(6)];
taug = t2(:);
