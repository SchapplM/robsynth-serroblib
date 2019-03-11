% Calculate potential energy for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:07
% EndTime: 2019-03-09 09:07:07
% DurationCPUTime: 0.72s
% Computational Cost: add. (253->95), mult. (483->102), div. (0->0), fcn. (530->10), ass. (0->48)
t71 = -m(6) - m(7);
t73 = pkin(9) - qJ(3);
t72 = -m(1) - m(2);
t29 = cos(pkin(6));
t33 = sin(qJ(1));
t36 = cos(qJ(2));
t56 = t33 * t36;
t32 = sin(qJ(2));
t37 = cos(qJ(1));
t57 = t32 * t37;
t14 = t29 * t57 + t56;
t28 = sin(pkin(6));
t52 = qJ(4) * t28;
t70 = t14 * pkin(3) + t37 * t52;
t69 = -mrSges(3,1) - mrSges(4,1) - mrSges(5,1);
t68 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t67 = -mrSges(5,3) + mrSges(4,2) + mrSges(3,3);
t66 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t65 = -m(7) * pkin(5) - t34 * mrSges(7,1) + t30 * mrSges(7,2) - mrSges(6,1);
t64 = t30 * mrSges(7,1) + t34 * mrSges(7,2) - t71 * t73 + mrSges(6,3) - t68;
t63 = t28 * t32;
t62 = t28 * t33;
t35 = cos(qJ(5));
t61 = t28 * t35;
t60 = t28 * t36;
t59 = t28 * t37;
t58 = t32 * t33;
t55 = t36 * t37;
t15 = t29 * t56 + t57;
t53 = qJ(3) * t15;
t51 = pkin(7) + r_base(3);
t50 = t29 * pkin(8) + t51;
t49 = t37 * pkin(1) + pkin(8) * t62 + r_base(1);
t48 = pkin(2) * t63 + t50;
t16 = -t29 * t58 + t55;
t47 = t16 * pkin(2) + t49;
t45 = t33 * pkin(1) - pkin(8) * t59 + r_base(2);
t44 = t14 * pkin(2) + t45;
t43 = pkin(3) * t63 - qJ(4) * t29 + t48;
t42 = t16 * pkin(3) - t33 * t52 + t47;
t13 = -t29 * t55 + t58;
t40 = t13 * qJ(3) + t44;
t38 = pkin(4) * t63 + t73 * t60 + t43;
t31 = sin(qJ(5));
t12 = -t29 * t31 + t32 * t61;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t50 - m(4) * t48 - m(5) * t43 - m(6) * t38 - t12 * mrSges(6,1) - mrSges(6,3) * t60 - m(7) * (pkin(5) * t12 + t38) - (t12 * t34 + t30 * t60) * mrSges(7,1) - (-t12 * t30 + t34 * t60) * mrSges(7,2) + t66 * (t29 * t35 + t31 * t63) - t67 * t29 + (((m(4) + m(5)) * qJ(3) + t68) * t36 + t69 * t32) * t28) * g(3) + (-m(3) * t45 - m(5) * (t40 + t70) - m(4) * t40 - mrSges(1,2) - t33 * mrSges(2,1) - t37 * mrSges(2,2) + t72 * r_base(2) + t71 * (t14 * pkin(4) + t44 + t70) + t66 * (t14 * t31 - t35 * t59) + t67 * t59 + t65 * (t14 * t35 + t31 * t59) + t69 * t14 + t64 * t13) * g(2) + (-m(4) * (t47 + t53) - m(5) * (t42 + t53) - m(3) * t49 - mrSges(1,1) + t33 * mrSges(2,2) - t37 * mrSges(2,1) + t72 * r_base(1) + t71 * (t16 * pkin(4) + t42) + t66 * (t16 * t31 + t33 * t61) - t67 * t62 + t65 * (t16 * t35 - t31 * t62) + t69 * t16 + t64 * t15) * g(1);
U  = t1;
