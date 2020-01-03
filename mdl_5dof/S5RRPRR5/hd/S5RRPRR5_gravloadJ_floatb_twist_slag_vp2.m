% Calculate Gravitation load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (275->64), mult. (207->68), div. (0->0), fcn. (157->10), ass. (0->40)
t38 = pkin(9) + qJ(4);
t30 = cos(t38);
t41 = cos(pkin(9));
t29 = sin(t38);
t31 = qJ(5) + t38;
t24 = sin(t31);
t57 = mrSges(6,2) * t24;
t48 = mrSges(5,2) * t29 + t57;
t68 = -mrSges(3,1) + sin(pkin(9)) * mrSges(4,2) + t48 - mrSges(5,1) * t30 - mrSges(4,1) * t41;
t67 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t52 = m(6) * pkin(4) + mrSges(5,1);
t66 = -mrSges(5,2) * t30 - t52 * t29;
t39 = qJ(1) + qJ(2);
t33 = cos(t39);
t25 = cos(t31);
t56 = t25 * t33;
t59 = mrSges(6,1) * t24;
t65 = mrSges(6,2) * t56 + t33 * t59;
t32 = sin(t39);
t64 = g(2) * t32;
t26 = t41 * pkin(3) + pkin(2);
t3 = pkin(4) * t30 + t26;
t42 = -pkin(7) - qJ(3);
t37 = -pkin(8) + t42;
t63 = t32 * t3 + t33 * t37;
t62 = t32 * t26 + t33 * t42;
t16 = t25 * mrSges(6,1);
t55 = t33 * pkin(2) + t32 * qJ(3);
t53 = -m(3) * pkin(1) - mrSges(2,1);
t51 = t33 * t3 - t32 * t37;
t50 = t33 * t26 - t32 * t42;
t49 = -mrSges(6,2) * t25 - t59;
t46 = -mrSges(6,1) * t56 + t67 * t32 + t68 * t33;
t45 = (m(4) * qJ(3) - t67) * t33 + (-t16 + t68) * t32;
t44 = cos(qJ(1));
t43 = sin(qJ(1));
t36 = t44 * pkin(1);
t35 = t43 * pkin(1);
t27 = t32 * pkin(2);
t1 = [(-t44 * mrSges(2,2) - m(4) * (t27 + t35) - m(5) * (t35 + t62) - m(6) * (t35 + t63) + t53 * t43 + t45) * g(3) + (t43 * mrSges(2,2) - m(4) * (t36 + t55) - m(5) * (t36 + t50) - m(6) * (t36 + t51) + t53 * t44 + t46) * g(2), (-m(4) * t27 - m(5) * t62 - m(6) * t63 + t45) * g(3) + (-m(4) * t55 - m(5) * t50 - m(6) * t51 + t46) * g(2), (g(2) * t33 + g(3) * t32) * (m(4) + m(5) + m(6)), (t66 * t33 - t65) * g(3) + (-t52 * t30 - t16 + t48) * g(1) + (-t49 - t66) * t64, -g(1) * (t16 - t57) - g(3) * t65 - t49 * t64];
taug = t1(:);
