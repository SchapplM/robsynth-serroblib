% Calculate potential energy for
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:24
% EndTime: 2019-03-09 12:27:25
% DurationCPUTime: 0.82s
% Computational Cost: add. (364->119), mult. (559->132), div. (0->0), fcn. (633->12), ass. (0->51)
t75 = -m(1) - m(2);
t74 = -mrSges(6,1) - mrSges(7,1);
t73 = -mrSges(6,2) - mrSges(7,2);
t72 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t71 = -mrSges(5,3) - t72;
t70 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t37 = sin(pkin(11));
t69 = pkin(3) * t37;
t40 = cos(pkin(6));
t47 = cos(qJ(2));
t48 = cos(qJ(1));
t59 = t47 * t48;
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t62 = t44 * t45;
t17 = -t40 * t59 + t62;
t43 = sin(qJ(5));
t68 = t17 * t43;
t60 = t45 * t47;
t61 = t44 * t48;
t19 = t40 * t60 + t61;
t67 = t19 * t43;
t38 = sin(pkin(6));
t66 = t38 * t44;
t65 = t38 * t45;
t64 = t38 * t47;
t63 = t38 * t48;
t58 = pkin(7) + r_base(3);
t57 = t45 * pkin(1) + r_base(2);
t56 = t37 * t65;
t55 = t43 * t64;
t54 = t40 * pkin(8) + t58;
t53 = t48 * pkin(1) + pkin(8) * t65 + r_base(1);
t52 = -pkin(8) * t63 + t57;
t39 = cos(pkin(11));
t29 = pkin(3) * t39 + pkin(2);
t42 = -pkin(9) - qJ(3);
t51 = t29 * t66 + t40 * t69 + t42 * t64 + t54;
t20 = -t40 * t62 + t59;
t50 = pkin(3) * t56 - t19 * t42 + t20 * t29 + t53;
t18 = t40 * t61 + t60;
t49 = t18 * t29 - t17 * t42 + (-pkin(8) - t69) * t63 + t57;
t46 = cos(qJ(5));
t36 = pkin(11) + qJ(4);
t32 = cos(t36);
t31 = sin(t36);
t30 = pkin(5) * t46 + pkin(4);
t14 = t31 * t40 + t32 * t66;
t10 = t20 * t32 + t31 * t65;
t8 = t18 * t32 - t31 * t63;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(5) * t51 - t14 * mrSges(5,1) + mrSges(5,3) * t64 - m(6) * (pkin(4) * t14 + t51) - m(7) * (-pkin(5) * t55 + t14 * t30 + t51) + t74 * (t14 * t46 - t55) + (-m(3) - m(4)) * t54 + t73 * (-t14 * t43 - t46 * t64) + (-t37 * mrSges(4,1) - t39 * mrSges(4,2) - mrSges(3,3)) * t40 + (t72 * t47 + (-m(4) * pkin(2) - t39 * mrSges(4,1) + t37 * mrSges(4,2) - mrSges(3,1)) * t44) * t38 + t70 * (t31 * t66 - t40 * t32)) * g(3) + (-m(4) * (t18 * pkin(2) + t52) - m(3) * t52 - m(6) * (t8 * pkin(4) + t49) - m(5) * t49 - m(7) * (pkin(5) * t68 + t8 * t30 + t49) + mrSges(3,3) * t63 - (-t18 * t37 - t39 * t63) * mrSges(4,2) - (t18 * t39 - t37 * t63) * mrSges(4,1) - mrSges(1,2) - t8 * mrSges(5,1) - t18 * mrSges(3,1) - t45 * mrSges(2,1) - t48 * mrSges(2,2) + t75 * r_base(2) + t74 * (t46 * t8 + t68) + t73 * (t17 * t46 - t43 * t8) + t70 * (t18 * t31 + t32 * t63) + t71 * t17) * g(2) + (-m(4) * (pkin(2) * t20 + t53) - m(3) * t53 - m(6) * (pkin(4) * t10 + t50) - m(5) * t50 - mrSges(3,3) * t65 - (-t20 * t37 + t39 * t65) * mrSges(4,2) - m(7) * (pkin(5) * t67 + t10 * t30 + t50) - (t20 * t39 + t56) * mrSges(4,1) - mrSges(1,1) - t10 * mrSges(5,1) - t20 * mrSges(3,1) + t45 * mrSges(2,2) - t48 * mrSges(2,1) + t75 * r_base(1) + t74 * (t10 * t46 + t67) + t73 * (-t10 * t43 + t19 * t46) + t70 * (t20 * t31 - t32 * t65) + t71 * t19) * g(1);
U  = t1;
