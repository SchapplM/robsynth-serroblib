% Calculate potential energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:19
% EndTime: 2019-12-05 16:02:20
% DurationCPUTime: 0.59s
% Computational Cost: add. (193->86), mult. (370->97), div. (0->0), fcn. (405->10), ass. (0->41)
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t58 = -mrSges(3,1) + mrSges(4,2);
t57 = mrSges(3,2) - mrSges(4,3);
t56 = mrSges(3,3) + mrSges(4,1);
t55 = m(6) * pkin(8) - mrSges(5,2) + mrSges(6,3);
t28 = sin(qJ(5));
t31 = cos(qJ(5));
t54 = -t28 * mrSges(6,1) - t31 * mrSges(6,2) - mrSges(5,3) + t58;
t53 = -m(6) * pkin(4) - mrSges(6,1) * t31 + mrSges(6,2) * t28 - mrSges(5,1);
t24 = sin(pkin(9));
t25 = sin(pkin(5));
t52 = t24 * t25;
t26 = cos(pkin(9));
t51 = t25 * t26;
t29 = sin(qJ(4));
t50 = t25 * t29;
t30 = sin(qJ(2));
t49 = t25 * t30;
t32 = cos(qJ(4));
t48 = t25 * t32;
t33 = cos(qJ(2));
t47 = t25 * t33;
t27 = cos(pkin(5));
t46 = t27 * t30;
t45 = t27 * t33;
t44 = pkin(6) * t51;
t43 = t24 * pkin(1) + r_base(2);
t42 = qJ(1) + r_base(3);
t41 = t26 * pkin(1) + pkin(6) * t52 + r_base(1);
t40 = t27 * pkin(6) + t42;
t39 = pkin(2) * t49 + t40;
t8 = t24 * t30 - t26 * t45;
t9 = t24 * t33 + t26 * t46;
t38 = t9 * pkin(2) + qJ(3) * t8 + t43;
t10 = t24 * t45 + t26 * t30;
t11 = -t24 * t46 + t26 * t33;
t37 = t11 * pkin(2) + qJ(3) * t10 + t41;
t36 = t27 * pkin(3) + pkin(7) * t49 - qJ(3) * t47 + t39;
t13 = t27 * t32 - t29 * t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t42 - mrSges(2,3) - m(3) * t40 - m(4) * t39 - m(5) * t36 - t13 * mrSges(5,1) - mrSges(5,3) * t49 - m(6) * (t13 * pkin(4) + t36) - (t13 * t31 + t28 * t49) * mrSges(6,1) - (-t13 * t28 + t31 * t49) * mrSges(6,2) - t56 * t27 + ((m(4) * qJ(3) - t57) * t33 + t58 * t30) * t25 - t55 * (t27 * t29 + t32 * t47)) * g(3) + (-mrSges(1,2) - t24 * mrSges(2,1) - t26 * mrSges(2,2) - m(3) * (t43 - t44) - m(4) * (t38 - t44) + t60 * r_base(2) + t57 * t8 + t56 * t51 + t59 * (pkin(7) * t9 + (-pkin(3) - pkin(6)) * t51 + t38) + t55 * (t26 * t50 + t32 * t8) + t53 * (-t26 * t48 + t29 * t8) + t54 * t9) * g(2) + (-m(3) * t41 - m(4) * t37 - t26 * mrSges(2,1) + t24 * mrSges(2,2) - mrSges(1,1) + t60 * r_base(1) - t56 * t52 + t59 * (pkin(3) * t52 + pkin(7) * t11 + t37) + t57 * t10 - t55 * (-t10 * t32 + t24 * t50) + t53 * (t10 * t29 + t24 * t48) + t54 * t11) * g(1);
U = t1;
