% Calculate Gravitation load on the joints for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:20
% EndTime: 2019-03-09 02:28:21
% DurationCPUTime: 0.40s
% Computational Cost: add. (218->83), mult. (309->97), div. (0->0), fcn. (250->8), ass. (0->47)
t27 = qJ(4) + qJ(5);
t22 = cos(t27);
t70 = (mrSges(6,2) - mrSges(7,3)) * t22;
t69 = -m(4) - m(5);
t65 = m(6) + m(7);
t68 = t65 * pkin(4) + mrSges(5,1);
t21 = sin(t27);
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t55 = t32 * mrSges(5,2);
t67 = t29 * mrSges(5,1) + t21 * mrSges(6,1) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t55 + t70;
t66 = m(5) * pkin(7) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t18 = t22 * pkin(9);
t64 = t29 * pkin(4);
t28 = sin(qJ(6));
t63 = mrSges(7,2) * t28;
t30 = sin(qJ(1));
t62 = t21 * t30;
t33 = cos(qJ(1));
t61 = t21 * t33;
t59 = t30 * t22;
t58 = t30 * t28;
t31 = cos(qJ(6));
t57 = t30 * t31;
t56 = t31 * mrSges(7,1);
t54 = t33 * t22;
t53 = t33 * t28;
t52 = t33 * t31;
t51 = -pkin(1) - qJ(3);
t49 = t33 * pkin(1) + t30 * qJ(2);
t48 = t22 * t56;
t47 = t33 * qJ(3) + t49;
t46 = -t65 + t69;
t44 = t51 - t64;
t42 = t21 * pkin(5) - t18;
t39 = mrSges(6,2) * t21 + t22 * t63;
t38 = -mrSges(6,1) * t59 - m(7) * (pkin(5) * t59 + pkin(9) * t62) - t30 * t48 - mrSges(7,3) * t62;
t37 = -mrSges(6,1) * t54 - m(7) * (pkin(5) * t54 + pkin(9) * t61) - t33 * t48 - mrSges(7,3) * t61;
t36 = -m(7) * t18 + (m(7) * pkin(5) + mrSges(6,1) + t56 - t63) * t21 + t70;
t35 = mrSges(5,2) * t29 - t68 * t32 + t39;
t34 = -pkin(8) - pkin(7);
t25 = t33 * qJ(2);
t4 = t21 * t52 - t58;
t3 = -t21 * t53 - t57;
t2 = -t21 * t57 - t53;
t1 = t21 * t58 - t52;
t5 = [(-m(3) * t49 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + t69 * t47 - t65 * (t30 * t34 + t33 * t64 + t47) + t66 * t30 + (-m(7) * t42 - t67) * t33) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) - t65 * (t33 * t34 + t25) + (-m(3) + t69) * t25 + t66 * t33 + (m(3) * pkin(1) - m(6) * t44 - m(7) * (-t42 + t44) + t69 * t51 + t67) * t30) * g(1) (-g(1) * t30 + g(2) * t33) * (m(3) - t46) (t33 * g(1) + t30 * g(2)) * t46 (t68 * t29 + t36 + t55) * g(3) + (t35 * t30 + t38) * g(2) + (t35 * t33 + t37) * g(1), t36 * g(3) + (t39 * t30 + t38) * g(2) + (t39 * t33 + t37) * g(1), -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t28 - mrSges(7,2) * t31) * t22];
taug  = t5(:);
