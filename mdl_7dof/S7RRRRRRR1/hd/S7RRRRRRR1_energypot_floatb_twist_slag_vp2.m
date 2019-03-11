% Calculate potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S7RRRRRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [8x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:06
% EndTime: 2019-03-10 06:29:07
% DurationCPUTime: 0.44s
% Computational Cost: add. (318->92), mult. (705->105), div. (0->0), fcn. (857->14), ass. (0->47)
t68 = -m(2) - m(3);
t67 = -m(4) - m(5);
t66 = -m(6) - m(7);
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t65 = t37 * t38;
t39 = sin(qJ(1));
t64 = t38 * t39;
t44 = cos(qJ(3));
t63 = t38 * t44;
t46 = cos(qJ(1));
t62 = t38 * t46;
t45 = cos(qJ(2));
t61 = t39 * t45;
t60 = t45 * t46;
t59 = mrSges(2,2) - mrSges(3,3);
t58 = -mrSges(3,2) - mrSges(4,3);
t57 = mrSges(4,2) + mrSges(5,3);
t56 = pkin(1) + r_base(3);
t55 = -m(1) - m(8) + t68;
t53 = -m(8) * pkin(3) + mrSges(5,2) - mrSges(6,3);
t52 = -m(8) * pkin(4) - mrSges(7,2) - mrSges(8,3);
t51 = -pkin(2) * t64 + r_base(2);
t50 = -pkin(2) * t62 + r_base(1);
t33 = sin(qJ(7));
t40 = cos(qJ(7));
t49 = -t40 * mrSges(8,1) + t33 * mrSges(8,2) - mrSges(7,1);
t48 = t33 * mrSges(8,1) + t40 * mrSges(8,2) + mrSges(6,2) - mrSges(7,3);
t47 = -t45 * mrSges(3,1) - mrSges(2,1) + (m(8) * pkin(2) - t58) * t38;
t43 = cos(qJ(4));
t42 = cos(qJ(5));
t41 = cos(qJ(6));
t36 = sin(qJ(4));
t35 = sin(qJ(5));
t34 = sin(qJ(6));
t25 = -t39 * t37 + t44 * t60;
t24 = -t37 * t60 - t39 * t44;
t23 = t37 * t46 + t44 * t61;
t22 = -t37 * t61 + t44 * t46;
t21 = -t36 * t45 + t43 * t63;
t18 = t25 * t43 + t36 * t62;
t17 = t25 * t36 - t43 * t62;
t16 = t23 * t43 + t36 * t64;
t15 = t23 * t36 - t43 * t64;
t10 = t18 * t42 + t24 * t35;
t8 = t16 * t42 + t22 * t35;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(1) * r_base(3) - t21 * mrSges(5,1) + t58 * t45 + t68 * t56 + t48 * (t21 * t35 + t42 * t65) + (-mrSges(4,1) * t44 + t57 * t37 - mrSges(3,1)) * t38 + (-t52 * t34 + t49 * t41 - mrSges(6,1)) * (t21 * t42 - t35 * t65) + (t66 * pkin(3) + t49 * t34 + t52 * t41 + t53) * (t36 * t63 + t43 * t45) + (-m(8) + t67 + t66) * (t45 * pkin(2) + t56)) * g(3) + (-t23 * mrSges(4,1) - t16 * mrSges(5,1) - t8 * mrSges(6,1) - mrSges(1,2) + t49 * (t15 * t34 + t41 * t8) + t48 * (t16 * t35 - t22 * t42) - t59 * t46 + t67 * t51 - t57 * t22 + t66 * (pkin(3) * t15 + t51) + t53 * t15 + t52 * (t15 * t41 - t34 * t8) + t55 * r_base(2) + t47 * t39) * g(2) + (-t25 * mrSges(4,1) - t18 * mrSges(5,1) - t10 * mrSges(6,1) - mrSges(1,1) + t49 * (t10 * t41 + t17 * t34) + t48 * (t18 * t35 - t24 * t42) + t59 * t39 + t67 * t50 - t57 * t24 + t52 * (-t10 * t34 + t17 * t41) + t53 * t17 + t66 * (t17 * pkin(3) + t50) + t55 * r_base(1) + t47 * t46) * g(1);
U  = t1;
