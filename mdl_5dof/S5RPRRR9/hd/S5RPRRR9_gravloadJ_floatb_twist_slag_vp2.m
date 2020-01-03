% Calculate Gravitation load on the joints for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:20
% DurationCPUTime: 0.37s
% Computational Cost: add. (283->76), mult. (271->86), div. (0->0), fcn. (226->10), ass. (0->42)
t27 = pkin(9) + qJ(3);
t24 = qJ(4) + t27;
t19 = sin(t24);
t20 = cos(t24);
t31 = sin(qJ(5));
t59 = t31 * mrSges(6,2);
t77 = -t19 * t59 + t20 * (-m(6) * pkin(8) - mrSges(6,3));
t71 = t20 * pkin(4) + t19 * pkin(8);
t76 = m(6) * t71;
t75 = -t20 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t19;
t72 = m(5) + m(6);
t33 = cos(qJ(5));
t56 = t33 * mrSges(6,1);
t70 = -(t56 - t59) * t20 + t75;
t30 = -pkin(6) - qJ(2);
t68 = -m(3) * qJ(2) + m(4) * t30 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t29 = cos(pkin(9));
t21 = t29 * pkin(2) + pkin(1);
t22 = sin(t27);
t23 = cos(t27);
t43 = t23 * mrSges(4,1) - t22 * mrSges(4,2);
t67 = m(4) * t21 + mrSges(2,1) + m(3) * pkin(1) + t29 * mrSges(3,1) - sin(pkin(9)) * mrSges(3,2) + t43 - t75;
t66 = pkin(3) * t22;
t18 = pkin(3) * t23;
t63 = mrSges(5,2) * t20;
t32 = sin(qJ(1));
t58 = t32 * t31;
t57 = t32 * t33;
t34 = cos(qJ(1));
t55 = t34 * t31;
t54 = t34 * t33;
t46 = t77 * t32;
t45 = t77 * t34;
t36 = m(6) * (-pkin(4) * t19 - t66) - t19 * t56;
t35 = t63 + (m(6) * pkin(4) + mrSges(5,1) + t56) * t19;
t26 = -pkin(7) + t30;
t8 = t18 + t21;
t4 = t20 * t54 + t58;
t3 = -t20 * t55 + t57;
t2 = -t20 * t57 + t55;
t1 = t20 * t58 + t54;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t72 * (-t32 * t26 + t34 * t8) + t68 * t32 + (-t67 - t76) * t34) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (t72 * t26 + t68) * t34 + (m(5) * t8 - m(6) * (-t71 - t8) + t67) * t32) * g(1), (-g(1) * t32 + g(2) * t34) * (m(3) + m(4) + t72), -g(1) * (t36 * t34 - t45) - g(2) * (t36 * t32 - t46) + (-t43 - m(5) * t18 - m(6) * (t18 + t71) + t70) * g(3) + (m(5) * t66 + mrSges(4,1) * t22 + mrSges(5,1) * t19 + mrSges(4,2) * t23 + t63) * (g(1) * t34 + g(2) * t32), (t70 - t76) * g(3) + (t35 * t32 + t46) * g(2) + (t35 * t34 + t45) * g(1), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t31 - mrSges(6,2) * t33) * t19];
taug = t5(:);
