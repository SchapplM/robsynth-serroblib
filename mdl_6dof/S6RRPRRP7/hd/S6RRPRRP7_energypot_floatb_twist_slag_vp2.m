% Calculate potential energy for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:39
% EndTime: 2019-03-09 12:16:40
% DurationCPUTime: 0.56s
% Computational Cost: add. (207->79), mult. (365->82), div. (0->0), fcn. (377->8), ass. (0->37)
t63 = -mrSges(3,1) - mrSges(4,1);
t62 = mrSges(3,2) - mrSges(4,3);
t61 = -m(1) - m(2);
t60 = -m(6) - m(7);
t32 = sin(qJ(2));
t36 = cos(qJ(2));
t59 = t62 * t32 + t63 * t36 - mrSges(2,1);
t58 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t57 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t56 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t55 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t35 = cos(qJ(4));
t54 = t32 * t35;
t33 = sin(qJ(1));
t53 = t33 * t36;
t37 = cos(qJ(1));
t52 = t36 * t37;
t51 = qJ(3) * t32;
t29 = pkin(6) + r_base(3);
t50 = t37 * pkin(1) + t33 * pkin(7) + r_base(1);
t31 = sin(qJ(4));
t12 = t31 * t32 + t35 * t36;
t47 = t33 * pkin(1) - pkin(7) * t37 + r_base(2);
t46 = pkin(2) * t52 + t37 * t51 + t50;
t45 = t32 * pkin(2) - qJ(3) * t36 + t29;
t44 = pkin(2) * t53 + t33 * t51 + t47;
t43 = t32 * pkin(3) + t45;
t42 = pkin(3) * t53 + t37 * pkin(8) + t44;
t41 = pkin(3) * t52 - pkin(8) * t33 + t46;
t34 = cos(qJ(5));
t30 = sin(qJ(5));
t13 = -t31 * t36 + t54;
t10 = t12 * t37;
t9 = t31 * t52 - t37 * t54;
t8 = t12 * t33;
t7 = t31 * t53 - t33 * t54;
t1 = (-m(1) * r_base(3) - m(4) * t45 - m(5) * t43 - mrSges(1,3) - mrSges(2,3) + t60 * (t13 * pkin(4) + pkin(9) * t12 + t43) - t62 * t36 + t63 * t32 + (-m(2) - m(3)) * t29 + (t55 * t30 + t56 * t34 - mrSges(5,1)) * t13 + t58 * t12) * g(3) + (-m(3) * t47 - m(4) * t44 - m(5) * t42 - t8 * mrSges(5,1) - mrSges(1,2) + t61 * r_base(2) + t60 * (t8 * pkin(4) + t7 * pkin(9) + t42) + t56 * (t30 * t37 + t8 * t34) + t55 * (t30 * t8 - t37 * t34) + t58 * t7 + t59 * t33 - t57 * t37) * g(2) + (-m(3) * t50 - m(4) * t46 - m(5) * t41 - t10 * mrSges(5,1) - mrSges(1,1) + t61 * r_base(1) + t60 * (t10 * pkin(4) + pkin(9) * t9 + t41) + t56 * (t10 * t34 - t30 * t33) + t55 * (t10 * t30 + t33 * t34) + t58 * t9 + t59 * t37 + t57 * t33) * g(1);
U  = t1;
