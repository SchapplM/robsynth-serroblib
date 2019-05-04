% Calculate potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10V2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:26
% EndTime: 2019-04-11 14:41:27
% DurationCPUTime: 0.56s
% Computational Cost: add. (244->79), mult. (291->87), div. (0->0), fcn. (295->12), ass. (0->39)
t60 = -m(2) - m(3);
t59 = -m(1) + t60;
t58 = -m(5) - m(6) - m(7);
t28 = qJ(2) + qJ(3);
t24 = sin(t28);
t25 = cos(t28);
t32 = sin(qJ(2));
t37 = cos(qJ(2));
t57 = -m(3) * pkin(1) - t37 * mrSges(3,1) - t25 * mrSges(4,1) + t32 * mrSges(3,2) + mrSges(4,2) * t24 - mrSges(2,1);
t56 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t55 = -m(7) * pkin(6) + mrSges(6,2) - mrSges(7,3);
t29 = sin(qJ(6));
t34 = cos(qJ(6));
t54 = -mrSges(7,1) * t34 + mrSges(7,2) * t29 - mrSges(6,1);
t53 = -mrSges(7,1) * t29 - mrSges(7,2) * t34 + mrSges(5,2) - mrSges(6,3);
t52 = pkin(3) * t25;
t31 = sin(qJ(4));
t51 = t24 * t31;
t36 = cos(qJ(4));
t50 = t24 * t36;
t33 = sin(qJ(1));
t49 = t33 * t24;
t48 = t33 * t31;
t47 = t33 * t36;
t38 = cos(qJ(1));
t46 = t38 * t24;
t45 = t38 * t31;
t44 = t38 * t36;
t27 = pkin(4) + r_base(3);
t23 = pkin(2) * t37 + pkin(1);
t43 = t23 * t33 + r_base(2);
t42 = t38 * t23 + r_base(1);
t41 = pkin(2) * t32 + t27;
t35 = cos(qJ(5));
t30 = sin(qJ(5));
t13 = t25 * t44 + t48;
t11 = t25 * t47 - t45;
t9 = -t25 * t30 + t35 * t50;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t32 * mrSges(3,1) - t37 * mrSges(3,2) - m(4) * t41 - t9 * mrSges(6,1) - mrSges(6,3) * t51 - (t29 * t51 + t9 * t34) * mrSges(7,1) - (-t9 * t29 + t34 * t51) * mrSges(7,2) + t55 * (t25 * t35 + t30 * t50) + t58 * (t24 * pkin(3) - t25 * pkin(5) + t41) + t60 * t27 + (-mrSges(4,2) + mrSges(5,3)) * t25 + (-mrSges(5,1) * t36 + mrSges(5,2) * t31 - mrSges(4,1)) * t24) * g(3) + (-m(4) * t43 - t11 * mrSges(5,1) - mrSges(5,3) * t49 - mrSges(1,2) + t54 * (t11 * t35 + t30 * t49) + t53 * (t25 * t48 + t44) + t58 * (pkin(5) * t49 + t33 * t52 + t43) + t55 * (t11 * t30 - t35 * t49) + t59 * r_base(2) - t56 * t38 + t57 * t33) * g(2) + (-m(4) * t42 - t13 * mrSges(5,1) - mrSges(5,3) * t46 - mrSges(1,1) + t58 * (pkin(5) * t46 + t38 * t52 + t42) + t55 * (t13 * t30 - t35 * t46) + t54 * (t13 * t35 + t30 * t46) + t53 * (t25 * t45 - t47) + t59 * r_base(1) + t57 * t38 + t56 * t33) * g(1);
U  = t1;
