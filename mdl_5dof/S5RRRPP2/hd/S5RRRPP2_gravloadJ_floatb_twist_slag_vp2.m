% Calculate Gravitation load on the joints for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:29
% EndTime: 2019-12-31 20:51:31
% DurationCPUTime: 0.45s
% Computational Cost: add. (260->63), mult. (270->64), div. (0->0), fcn. (215->6), ass. (0->31)
t79 = mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t78 = -mrSges(4,1) - mrSges(5,1);
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t76 = t79 * t26 - t78 * t28;
t25 = qJ(1) + qJ(2);
t20 = sin(t25);
t21 = cos(t25);
t65 = g(1) * t21 + g(2) * t20;
t74 = -t28 * mrSges(6,1) - t76;
t73 = mrSges(6,3) - mrSges(5,2) - mrSges(4,3) + mrSges(3,2);
t67 = m(5) + m(6);
t23 = t28 * pkin(3);
t41 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t22 = t26 * qJ(4);
t45 = -pkin(2) - t22;
t63 = t73 * t21 + (-m(6) * t45 - t41 * t28 - m(5) * (t45 - t23) + mrSges(3,1) + t76) * t20;
t62 = (-mrSges(3,1) + t74) * t21 + t73 * t20;
t27 = sin(qJ(1));
t59 = t27 * pkin(1);
t29 = cos(qJ(1));
t24 = t29 * pkin(1);
t58 = t21 * t28;
t51 = t21 * pkin(2) + t20 * pkin(7);
t50 = t23 + t22;
t18 = t21 * pkin(7);
t47 = -t20 * pkin(2) + t18;
t44 = pkin(3) * t58 + t21 * t22 + t51;
t43 = -t21 * qJ(5) + t18;
t34 = pkin(4) * t58 - t20 * qJ(5) + t44;
t1 = [(-t29 * mrSges(2,1) + t27 * mrSges(2,2) - m(3) * t24 - m(4) * (t24 + t51) - m(5) * (t24 + t44) - m(6) * (t24 + t34) + t62) * g(2) + (t27 * mrSges(2,1) + t29 * mrSges(2,2) + m(3) * t59 - m(4) * (t47 - t59) - m(5) * (t18 - t59) - m(6) * (t43 - t59) + t63) * g(1), (-m(4) * t51 - m(5) * t44 - m(6) * t34 + t62) * g(2) + (-m(4) * t47 - m(5) * t18 - m(6) * t43 + t63) * g(1), (-m(5) * t50 - m(6) * (t28 * pkin(4) + t50) + t74) * g(3) + ((m(5) * pkin(3) - t41 - t78) * t26 + (-qJ(4) * t67 - t79) * t28) * t65, (g(3) * t28 - t26 * t65) * t67, (g(1) * t20 - g(2) * t21) * m(6)];
taug = t1(:);
