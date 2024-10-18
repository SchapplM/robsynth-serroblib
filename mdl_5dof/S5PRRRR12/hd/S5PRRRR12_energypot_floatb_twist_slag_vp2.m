% Calculate potential energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:08
% EndTime: 2024-09-28 18:07:09
% DurationCPUTime: 0.41s
% Computational Cost: add. (334->133), mult. (410->165), div. (0->0), fcn. (399->26), ass. (0->69)
t41 = qJ(2) + qJ(3);
t38 = qJ(4) + t41;
t29 = cos(t38);
t46 = cos(pkin(6));
t93 = t29 * t46;
t32 = pkin(5) + t41;
t26 = qJ(4) + t32;
t16 = sin(t26) / 0.2e1;
t33 = pkin(5) - t41;
t27 = -qJ(4) + t33;
t17 = cos(t27) / 0.2e1;
t19 = sin(t32) / 0.2e1;
t20 = cos(t33) / 0.2e1;
t22 = sin(t27);
t23 = cos(t26);
t24 = sin(t33);
t25 = cos(t32);
t28 = sin(t38);
t56 = pkin(7) + pkin(8);
t40 = pkin(9) + t56;
t44 = sin(pkin(5));
t47 = cos(pkin(5));
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t55 = cos(qJ(2));
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t54 = cos(qJ(3));
t59 = pkin(3) * t50 * t55 + (pkin(3) * t54 + pkin(2)) * t51;
t70 = t47 * t52;
t71 = t47 * t51;
t72 = t47 * t48;
t43 = sin(pkin(6));
t79 = t43 * t44;
t85 = m(3) * pkin(7) + m(6) * (pkin(10) * t46 + t40) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t92 = (t70 * t28 - t48 * t79 + t72 * t93) * mrSges(6,1) - (t72 * t28 + t52 * t79 - t70 * t93) * mrSges(6,2) + m(4) * (pkin(2) * t71 - t44 * t56) + m(5) * (-t40 * t44 + t47 * t59) + t71 * mrSges(3,1) + (t19 - t24 / 0.2e1) * mrSges(4,1) + (t16 - t22 / 0.2e1) * mrSges(5,1) + t47 * t55 * mrSges(3,2) + (t20 + t25 / 0.2e1) * mrSges(4,2) + (t17 + t23 / 0.2e1) * mrSges(5,2) - t85 * t44 + mrSges(2,2);
t31 = t55 * pkin(2) + pkin(1);
t35 = cos(t41);
t73 = t46 * t52;
t74 = t46 * t48;
t88 = -mrSges(2,1) + t51 * mrSges(3,2) - m(4) * t31 - t35 * mrSges(4,1) + sin(t41) * mrSges(4,2) - m(5) * (pkin(3) * t35 + t31) - t29 * mrSges(5,1) - (-t74 * t28 + t52 * t29) * mrSges(6,1) - (-t73 * t28 - t48 * t29) * mrSges(6,2) + (-m(6) * pkin(2) - mrSges(3,1)) * t55;
t87 = -m(2) - m(4) - m(5);
t86 = -m(1) + t87;
t84 = pkin(10) * t43;
t42 = sin(pkin(11));
t83 = t28 * t42;
t82 = t29 * t43;
t80 = t42 * t47;
t45 = cos(pkin(11));
t78 = t43 * t45;
t76 = t44 * t46;
t75 = t45 * t47;
t68 = t42 * t84;
t67 = pkin(10) * t78;
t66 = t42 * pkin(1) + r_base(2);
t65 = t45 * pkin(1) + r_base(1);
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t7 = -pkin(4) * t42 + t47 * t67;
t9 = pkin(4) * t75 + t68;
t58 = pkin(3) * t42 + t49 * t9 - t53 * t7;
t6 = pkin(4) * t45 + t47 * t68;
t8 = pkin(4) * t80 - t67;
t57 = pkin(3) * t45 - t49 * t8 + t53 * t6;
t14 = -pkin(4) * t49 + t53 * t84;
t13 = pkin(4) * t53 + t49 * t84 + pkin(3);
t2 = pkin(3) * t75 + t49 * t7 + t53 * t9;
t1 = -pkin(3) * t80 - t49 * t6 - t53 * t8;
t3 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - (t20 - t25 / 0.2e1) * mrSges(4,1) - (t19 + t24 / 0.2e1) * mrSges(4,2) - (t17 - t23 / 0.2e1) * mrSges(5,1) - (t16 + t22 / 0.2e1) * mrSges(5,2) + (-m(5) * t59 - (t28 * t52 + t29 * t74) * mrSges(6,1) - (-t28 * t48 + t29 * t73) * mrSges(6,2) + mrSges(6,3) * t82 + (-mrSges(3,2) + m(6) * (-t13 * t50 + t14 * t54)) * t55 + (-mrSges(3,1) - m(4) * pkin(2) - m(6) * (t13 * t54 + t14 * t50 + pkin(2))) * t51) * t44 + (-m(3) - m(6) + t87) * (qJ(1) + r_base(3)) + (-m(4) * t56 - m(5) * t40 - t46 * mrSges(6,3) - (mrSges(6,1) * t48 + mrSges(6,2) * t52) * t43 - t85) * t47) * g(3) + (-mrSges(1,2) - m(3) * t66 + t83 * mrSges(5,2) - m(6) * ((t2 * t50 + t54 * t58) * t55 + (pkin(2) * t75 + t2 * t54 - t50 * t58) * t51 + t66) - (-t29 * t75 + t83) * t43 * mrSges(6,3) + t86 * r_base(2) + t88 * t42 + (t76 * mrSges(6,3) - t92) * t45) * g(2) + (-mrSges(1,1) - m(3) * t65 - m(6) * ((t1 * t50 + t54 * t57) * t55 + (-pkin(2) * t80 + t1 * t54 - t50 * t57) * t51 + t65) - t28 * t78 * mrSges(6,3) + t86 * r_base(1) + (t28 * mrSges(5,2) + t88) * t45 + (-(t47 * t82 + t76) * mrSges(6,3) + t92) * t42) * g(1);
U = t3;
