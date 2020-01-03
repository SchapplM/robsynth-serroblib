% Calculate Gravitation load on the joints for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (172->65), mult. (231->78), div. (0->0), fcn. (198->8), ass. (0->38)
t21 = qJ(2) + qJ(3);
t18 = sin(t21);
t19 = cos(t21);
t22 = sin(qJ(4));
t50 = t22 * mrSges(5,2);
t68 = -t18 * t50 + t19 * (-m(5) * pkin(7) - mrSges(5,3));
t62 = t19 * pkin(3) + t18 * pkin(7);
t67 = m(5) * t62;
t66 = -t19 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t18;
t63 = m(4) + m(5);
t25 = cos(qJ(4));
t46 = t25 * mrSges(5,1);
t61 = -(t46 - t50) * t19 + t66;
t59 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t36 = t26 * mrSges(3,1) - t23 * mrSges(3,2);
t58 = m(3) * pkin(1) + mrSges(2,1) + t36 - t66;
t57 = pkin(2) * t23;
t20 = t26 * pkin(2);
t54 = mrSges(4,2) * t19;
t24 = sin(qJ(1));
t49 = t22 * t24;
t27 = cos(qJ(1));
t48 = t22 * t27;
t47 = t24 * t25;
t45 = t25 * t27;
t39 = t68 * t24;
t38 = t68 * t27;
t30 = m(5) * (-pkin(3) * t18 - t57) - t18 * t46;
t29 = t54 + (m(5) * pkin(3) + mrSges(4,1) + t46) * t18;
t28 = -pkin(6) - pkin(5);
t17 = t20 + pkin(1);
t4 = t19 * t45 + t49;
t3 = -t19 * t48 + t47;
t2 = -t19 * t47 + t48;
t1 = t19 * t49 + t45;
t5 = [(-t4 * mrSges(5,1) - t3 * mrSges(5,2) - t63 * (t27 * t17 - t24 * t28) + t59 * t24 + (-t58 - t67) * t27) * g(2) + (-t2 * mrSges(5,1) - t1 * mrSges(5,2) + (t63 * t28 + t59) * t27 + (m(4) * t17 - m(5) * (-t17 - t62) + t58) * t24) * g(1), -g(1) * (t30 * t27 - t38) - g(2) * (t30 * t24 - t39) + (-t36 - m(4) * t20 - m(5) * (t20 + t62) + t61) * g(3) + (m(4) * t57 + mrSges(3,1) * t23 + mrSges(4,1) * t18 + mrSges(3,2) * t26 + t54) * (g(1) * t27 + g(2) * t24), (t61 - t67) * g(3) + (t29 * t24 + t39) * g(2) + (t29 * t27 + t38) * g(1), -g(1) * (mrSges(5,1) * t3 - mrSges(5,2) * t4) - g(2) * (-mrSges(5,1) * t1 + mrSges(5,2) * t2) - g(3) * (-mrSges(5,1) * t22 - mrSges(5,2) * t25) * t18];
taug = t5(:);
