% Calculate potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:01
% EndTime: 2019-03-09 10:38:02
% DurationCPUTime: 0.71s
% Computational Cost: add. (393->103), mult. (832->115), div. (0->0), fcn. (1011->12), ass. (0->59)
t78 = -mrSges(3,3) - mrSges(4,3);
t35 = sin(pkin(11));
t40 = sin(qJ(2));
t44 = cos(qJ(2));
t59 = cos(pkin(11));
t24 = -t40 * t35 + t44 * t59;
t77 = -m(3) - m(2) - m(1);
t76 = -m(3) * pkin(1) - mrSges(2,1);
t75 = m(3) * pkin(8) - t78;
t74 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t73 = -t38 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t72 = m(7) * (pkin(5) + pkin(9)) + t42 * mrSges(7,1) - t38 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t23 = -t44 * t35 - t40 * t59;
t41 = sin(qJ(1));
t45 = cos(qJ(1));
t37 = cos(pkin(6));
t48 = t37 * t24;
t9 = t41 * t23 + t45 * t48;
t71 = pkin(9) * t9;
t69 = pkin(2) * t40;
t11 = t23 * t45 - t41 * t48;
t68 = pkin(9) * t11;
t36 = sin(pkin(6));
t19 = t24 * t36;
t67 = pkin(9) * t19;
t66 = t36 * t41;
t65 = t36 * t45;
t63 = t40 * t45;
t62 = t41 * t40;
t61 = t41 * t44;
t60 = t44 * t45;
t58 = pkin(7) + r_base(3);
t57 = t37 * pkin(8) + t58;
t22 = t37 * t69 + (-pkin(8) - qJ(3)) * t36;
t32 = pkin(2) * t44 + pkin(1);
t55 = t45 * t22 + t41 * t32 + r_base(2);
t21 = t23 * t37;
t10 = -t21 * t45 + t41 * t24;
t54 = t10 * pkin(3) + t55;
t53 = -t22 * t41 + t45 * t32 + r_base(1);
t52 = t37 * qJ(3) + t36 * t69 + t57;
t12 = t41 * t21 + t24 * t45;
t51 = t12 * pkin(3) + t53;
t20 = t23 * t36;
t50 = -t20 * pkin(3) + t52;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t3 = t10 * t39 + t43 * t65;
t4 = t10 * t43 - t39 * t65;
t49 = t4 * pkin(4) + qJ(5) * t3 + t54;
t5 = t12 * t39 - t43 * t66;
t6 = t12 * t43 + t39 * t66;
t47 = t6 * pkin(4) + qJ(5) * t5 + t51;
t14 = -t20 * t39 - t37 * t43;
t15 = -t20 * t43 + t37 * t39;
t46 = t15 * pkin(4) + qJ(5) * t14 + t50;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t57 - (t40 * mrSges(3,1) + t44 * mrSges(3,2)) * t36 - m(4) * t52 + t20 * mrSges(4,1) - m(5) * (t50 - t67) - m(6) * (t46 - t67) - m(7) * t46 + t78 * t37 + t73 * t14 + t72 * t19 + t74 * t15) * g(3) + (-m(5) * (t54 - t71) - m(6) * (t49 - t71) - m(7) * t49 - (t37 * t63 + t61) * mrSges(3,1) - (t37 * t60 - t62) * mrSges(3,2) - m(4) * t55 - mrSges(1,2) - t10 * mrSges(4,1) - t45 * mrSges(2,2) + t76 * t41 + t77 * r_base(2) + t75 * t65 + t73 * t3 + t72 * t9 + t74 * t4) * g(2) + (-m(5) * (t51 - t68) - m(6) * (t47 - t68) - m(7) * t47 - (-t37 * t61 - t63) * mrSges(3,2) - (-t37 * t62 + t60) * mrSges(3,1) - m(4) * t53 - mrSges(1,1) - t12 * mrSges(4,1) + t41 * mrSges(2,2) + t76 * t45 + t77 * r_base(1) - t75 * t66 + t74 * t6 + t73 * t5 + t72 * t11) * g(1);
U  = t1;
