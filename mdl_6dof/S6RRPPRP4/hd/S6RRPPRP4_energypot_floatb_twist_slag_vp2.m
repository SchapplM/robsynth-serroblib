% Calculate potential energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:55
% EndTime: 2019-03-09 08:36:56
% DurationCPUTime: 0.61s
% Computational Cost: add. (224->85), mult. (406->91), div. (0->0), fcn. (430->8), ass. (0->38)
t67 = -m(1) - m(2);
t66 = -m(6) - m(7);
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t65 = -t38 * mrSges(3,1) + t35 * mrSges(3,2) - mrSges(2,1);
t64 = -mrSges(4,1) - mrSges(5,1);
t63 = mrSges(2,2) - mrSges(3,3);
t62 = mrSges(4,2) - mrSges(5,3);
t61 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3) + mrSges(7,2);
t60 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t59 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t32 = sin(pkin(9));
t58 = t32 * t35;
t33 = cos(pkin(9));
t57 = t33 * t35;
t36 = sin(qJ(1));
t56 = t35 * t36;
t39 = cos(qJ(1));
t55 = t35 * t39;
t54 = t36 * t38;
t53 = t38 * t39;
t52 = qJ(3) * t35;
t31 = pkin(6) + r_base(3);
t51 = t39 * pkin(1) + t36 * pkin(7) + r_base(1);
t49 = t36 * pkin(1) - pkin(7) * t39 + r_base(2);
t48 = pkin(2) * t53 + t39 * t52 + t51;
t47 = t35 * pkin(2) - qJ(3) * t38 + t31;
t46 = pkin(2) * t54 + t36 * t52 + t49;
t45 = pkin(3) * t57 + qJ(4) * t58 + t47;
t15 = t32 * t53 - t36 * t33;
t16 = t36 * t32 + t33 * t53;
t44 = t16 * pkin(3) + t15 * qJ(4) + t48;
t13 = t32 * t54 + t33 * t39;
t14 = -t32 * t39 + t33 * t54;
t42 = t14 * pkin(3) + t13 * qJ(4) + t46;
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - m(4) * t47 - m(5) * t45 - mrSges(1,3) - mrSges(2,3) + t66 * (pkin(4) * t57 + t38 * pkin(8) + t45) + t59 * (t34 * t57 - t37 * t58) + (-m(2) - m(3)) * t31 + (-mrSges(3,2) - t61) * t38 + (t60 * (t32 * t34 + t33 * t37) + t62 * t32 + t64 * t33 - mrSges(3,1)) * t35) * g(3) + (-m(3) * t49 - m(4) * t46 - m(5) * t42 - mrSges(1,2) + t67 * r_base(2) + t66 * (t14 * pkin(4) - pkin(8) * t56 + t42) - t63 * t39 + t65 * t36 + t60 * (t13 * t34 + t14 * t37) + t64 * t14 + t62 * t13 + t59 * (-t13 * t37 + t14 * t34) + t61 * t56) * g(2) + (-m(3) * t51 - m(4) * t48 - m(5) * t44 - mrSges(1,1) + t67 * r_base(1) + t66 * (t16 * pkin(4) - pkin(8) * t55 + t44) + t60 * (t15 * t34 + t16 * t37) + t65 * t39 + t63 * t36 + t59 * (-t15 * t37 + t16 * t34) + t64 * t16 + t62 * t15 + t61 * t55) * g(1);
U  = t1;
