% Calculate Gravitation load on the joints for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:42
% EndTime: 2019-12-31 20:05:44
% DurationCPUTime: 0.64s
% Computational Cost: add. (273->78), mult. (412->88), div. (0->0), fcn. (382->8), ass. (0->44)
t69 = mrSges(5,1) + mrSges(6,1);
t68 = -mrSges(5,2) + mrSges(6,3);
t67 = m(5) + m(6);
t22 = sin(qJ(2));
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t63 = g(1) * t25 + g(2) * t23;
t66 = t63 * t22;
t65 = -mrSges(5,3) - mrSges(6,2);
t18 = pkin(8) + qJ(4);
t13 = sin(t18);
t14 = cos(t18);
t64 = t68 * t13 + t69 * t14;
t62 = -m(3) - m(4);
t24 = cos(qJ(2));
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t30 = m(4) * pkin(2) + t20 * mrSges(4,1) - t19 * mrSges(4,2);
t60 = t30 * t24;
t59 = -t19 * mrSges(4,1) - t20 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t37 = t24 * mrSges(3,1) - t22 * mrSges(3,2);
t58 = t65 * t22 - t37;
t56 = m(6) * pkin(4) + t69;
t55 = m(6) * qJ(5) + t68;
t54 = pkin(3) * t19;
t51 = g(3) * t22;
t21 = -pkin(7) - qJ(3);
t48 = t22 * t21;
t47 = t23 * t24;
t12 = t20 * pkin(3) + pkin(2);
t6 = t24 * t12;
t46 = t25 * t13;
t45 = t25 * t14;
t44 = t25 * t22;
t43 = t25 * pkin(1) + t23 * pkin(6);
t38 = m(4) * qJ(3) + mrSges(4,3);
t32 = pkin(4) * t14 + qJ(5) * t13;
t29 = t38 * t22 + t60;
t16 = t25 * pkin(6);
t4 = t23 * t13 + t24 * t45;
t3 = -t23 * t14 + t24 * t46;
t2 = t14 * t47 - t46;
t1 = t13 * t47 + t45;
t5 = [(t65 * t44 + t62 * t43 - t67 * (-t21 * t44 + t23 * t54 + t25 * t6 + t43) - t56 * t4 - t55 * t3 + (-mrSges(2,1) - t37 - t29) * t25 + t59 * t23) * g(2) + (-t67 * (t23 * t48 + t25 * t54 + t16) + t56 * t2 + t62 * t16 + t55 * t1 + t59 * t25 + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (-t22 * qJ(3) - pkin(1)) + t22 * mrSges(4,3) + t60 - t67 * (-pkin(1) - t6) - t58) * t23) * g(1), (-t29 - t67 * (t6 - t48) + t58) * g(3) + ((-m(6) * t32 - t64) * g(3) + t63 * (t67 * t21 + mrSges(3,2) - t38 + t65)) * t24 + (mrSges(3,1) + t30 + m(5) * t12 - m(6) * (-t12 - t32) + t64) * t66, (t24 * g(3) - t66) * (m(4) + t67), (t56 * t13 - t55 * t14) * t51 + (t56 * t1 - t55 * t2) * g(2) + (t56 * t3 - t55 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - t13 * t51) * m(6)];
taug = t5(:);
