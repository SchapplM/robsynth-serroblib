% Calculate Gravitation load on the joints for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:13
% EndTime: 2019-03-09 03:11:15
% DurationCPUTime: 0.68s
% Computational Cost: add. (351->76), mult. (418->85), div. (0->0), fcn. (375->8), ass. (0->39)
t67 = -m(6) - m(7);
t60 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t59 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t75 = mrSges(4,1) - mrSges(5,2);
t74 = -mrSges(4,2) + mrSges(5,3);
t73 = -mrSges(7,2) - mrSges(6,3);
t23 = sin(qJ(5));
t26 = cos(qJ(5));
t72 = -t60 * t23 + t59 * t26;
t71 = -t73 + t67 * (-pkin(3) - pkin(8));
t22 = qJ(1) + pkin(9);
t16 = sin(t22);
t17 = cos(t22);
t65 = g(1) * t17 + g(2) * t16;
t24 = sin(qJ(3));
t27 = cos(qJ(3));
t66 = t74 * t24 + t75 * t27;
t46 = m(5) - t67;
t64 = mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t62 = t73 * t27 - t66;
t25 = sin(qJ(1));
t57 = pkin(1) * t25;
t54 = g(3) * t27;
t20 = t27 * pkin(3);
t28 = cos(qJ(1));
t21 = t28 * pkin(1);
t53 = t17 * t27;
t52 = t23 * t24;
t51 = t24 * t26;
t18 = t24 * qJ(4);
t48 = t20 + t18;
t45 = t17 * pkin(2) + t16 * pkin(7) + t21;
t43 = t17 * pkin(7) - t57;
t38 = pkin(3) * t53 + t17 * t18 + t45;
t4 = -t16 * t52 + t17 * t26;
t3 = t16 * t51 + t17 * t23;
t2 = t16 * t26 + t17 * t52;
t1 = t16 * t23 - t17 * t51;
t5 = [(-m(3) * t21 - m(4) * t45 - m(5) * t38 - t28 * mrSges(2,1) + t25 * mrSges(2,2) + t67 * (t16 * pkin(4) + pkin(8) * t53 + t38) - t60 * t2 - t59 * t1 + t64 * t16 + (-mrSges(3,1) + t62) * t17) * g(2) + (m(3) * t57 + t25 * mrSges(2,1) + t28 * mrSges(2,2) + (-m(4) - m(5)) * t43 + t67 * (t17 * pkin(4) + t43) - t60 * t4 - t59 * t3 + t64 * t17 + (m(4) * pkin(2) + m(5) * t20 + mrSges(3,1) - t46 * (-pkin(2) - t18) + t71 * t27 + t66) * t16) * g(1) (-m(3) - m(4) - t46) * g(3) (-m(5) * t48 + t67 * (t27 * pkin(8) + t48) + t72 * t24 + t62) * g(3) + ((m(5) * pkin(3) + t71 + t75) * t24 + (-qJ(4) * t46 + t72 - t74) * t27) * t65 (-t65 * t24 + t54) * t46 (t59 * t23 + t60 * t26) * t54 + (-t60 * t3 + t59 * t4) * g(2) + (t60 * t1 - t59 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t26 * t54) * m(7)];
taug  = t5(:);
