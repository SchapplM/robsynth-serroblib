% Calculate potential energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:03
% EndTime: 2019-03-09 09:17:04
% DurationCPUTime: 0.70s
% Computational Cost: add. (253->92), mult. (483->101), div. (0->0), fcn. (530->10), ass. (0->44)
t68 = -m(1) - m(2);
t67 = -m(6) - m(7);
t66 = -mrSges(4,3) - mrSges(5,1) + mrSges(3,2);
t65 = -mrSges(5,3) + mrSges(4,2) + mrSges(3,3);
t64 = -mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t63 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t62 = -m(7) * pkin(5) - t34 * mrSges(7,1) + t30 * mrSges(7,2) - mrSges(6,1);
t61 = -t30 * mrSges(7,1) - t34 * mrSges(7,2) + t64;
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t60 = t28 * t32;
t33 = sin(qJ(1));
t59 = t28 * t33;
t36 = cos(qJ(2));
t58 = t28 * t36;
t37 = cos(qJ(1));
t57 = t28 * t37;
t56 = t32 * t33;
t55 = t32 * t37;
t54 = t33 * t36;
t53 = t36 * t37;
t52 = qJ(4) * t28;
t51 = pkin(7) + r_base(3);
t29 = cos(pkin(6));
t49 = t29 * pkin(8) + t51;
t48 = t37 * pkin(1) + pkin(8) * t59 + r_base(1);
t47 = pkin(2) * t60 + t49;
t46 = t33 * pkin(1) - pkin(8) * t57 + r_base(2);
t15 = t29 * t54 + t55;
t16 = -t29 * t56 + t53;
t45 = t16 * pkin(2) + qJ(3) * t15 + t48;
t44 = pkin(3) * t60 - qJ(4) * t29 + t47;
t43 = pkin(9) * t60 + t44;
t13 = -t29 * t53 + t56;
t14 = t29 * t55 + t54;
t42 = t14 * pkin(2) + t13 * qJ(3) + t46;
t41 = t14 * pkin(3) + t37 * t52 + t42;
t40 = t16 * pkin(3) - t33 * t52 + t45;
t35 = cos(qJ(5));
t31 = sin(qJ(5));
t12 = -t29 * t31 - t35 * t58;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t49 - m(4) * t47 - m(5) * t44 - m(6) * t43 - t12 * mrSges(6,1) - m(7) * (pkin(5) * t12 + t43) - (t12 * t34 + t30 * t60) * mrSges(7,1) - (-t12 * t30 + t34 * t60) * mrSges(7,2) + t63 * (-t29 * t35 + t31 * t58) - t65 * t29 + (t64 * t32 + (t67 * (-pkin(4) - qJ(3)) + (m(4) + m(5)) * qJ(3) - t66) * t36) * t28) * g(3) + (-m(3) * t46 - m(4) * t42 - m(5) * t41 - t33 * mrSges(2,1) - t37 * mrSges(2,2) - mrSges(1,2) + t68 * r_base(2) + t67 * (t13 * pkin(4) + t14 * pkin(9) + t41) - t63 * (t13 * t31 - t35 * t57) + t65 * t57 + t62 * (t13 * t35 + t31 * t57) + t66 * t13 + t61 * t14) * g(2) + (-m(3) * t48 - m(4) * t45 - m(5) * t40 - t37 * mrSges(2,1) + t33 * mrSges(2,2) - mrSges(1,1) + t68 * r_base(1) + t67 * (t15 * pkin(4) + pkin(9) * t16 + t40) - t63 * (t15 * t31 + t35 * t59) - t65 * t59 + t62 * (t15 * t35 - t31 * t59) + t66 * t15 + t61 * t16) * g(1);
U  = t1;
