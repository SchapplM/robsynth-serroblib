% Calculate potential energy for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR16_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:37
% EndTime: 2019-12-31 20:44:37
% DurationCPUTime: 0.57s
% Computational Cost: add. (193->86), mult. (370->93), div. (0->0), fcn. (405->10), ass. (0->41)
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t58 = -mrSges(3,1) + mrSges(4,2);
t57 = mrSges(3,2) - mrSges(4,3);
t56 = mrSges(3,3) + mrSges(4,1);
t55 = m(6) * pkin(9) - mrSges(5,2) + mrSges(6,3);
t26 = sin(qJ(5));
t30 = cos(qJ(5));
t54 = -mrSges(6,1) * t26 - mrSges(6,2) * t30 - mrSges(5,3) + t58;
t53 = -m(6) * pkin(4) - mrSges(6,1) * t30 + mrSges(6,2) * t26 - mrSges(5,1);
t24 = sin(pkin(5));
t28 = sin(qJ(2));
t52 = t24 * t28;
t29 = sin(qJ(1));
t51 = t24 * t29;
t32 = cos(qJ(2));
t50 = t24 * t32;
t33 = cos(qJ(1));
t49 = t24 * t33;
t48 = t28 * t29;
t47 = t28 * t33;
t46 = t29 * t32;
t45 = t32 * t33;
t44 = pkin(6) + r_base(3);
t43 = pkin(7) * t49;
t42 = t29 * pkin(1) + r_base(2);
t25 = cos(pkin(5));
t41 = t25 * pkin(7) + t44;
t40 = t33 * pkin(1) + pkin(7) * t51 + r_base(1);
t39 = pkin(2) * t52 + t41;
t10 = -t25 * t45 + t48;
t11 = t25 * t47 + t46;
t38 = t11 * pkin(2) + t10 * qJ(3) + t42;
t12 = t25 * t46 + t47;
t13 = -t25 * t48 + t45;
t37 = t13 * pkin(2) + qJ(3) * t12 + t40;
t36 = t25 * pkin(3) + pkin(8) * t52 - qJ(3) * t50 + t39;
t31 = cos(qJ(4));
t27 = sin(qJ(4));
t9 = t25 * t31 - t27 * t50;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t44 - mrSges(2,3) - m(3) * t41 - m(4) * t39 - m(5) * t36 - t9 * mrSges(5,1) - mrSges(5,3) * t52 - m(6) * (pkin(4) * t9 + t36) - (t26 * t52 + t30 * t9) * mrSges(6,1) - (-t26 * t9 + t30 * t52) * mrSges(6,2) - t55 * (t25 * t27 + t31 * t50) - t56 * t25 + ((m(4) * qJ(3) - t57) * t32 + t58 * t28) * t24) * g(3) + (-mrSges(1,2) - t29 * mrSges(2,1) - t33 * mrSges(2,2) - m(3) * (t42 - t43) - m(4) * (t38 - t43) + t60 * r_base(2) + t56 * t49 + t59 * (t11 * pkin(8) + (-pkin(3) - pkin(7)) * t49 + t38) + t55 * (t10 * t31 + t27 * t49) + t57 * t10 + t53 * (t10 * t27 - t31 * t49) + t54 * t11) * g(2) + (-m(3) * t40 - m(4) * t37 - t33 * mrSges(2,1) + t29 * mrSges(2,2) - mrSges(1,1) + t60 * r_base(1) - t56 * t51 + t59 * (pkin(3) * t51 + pkin(8) * t13 + t37) + t57 * t12 - t55 * (-t12 * t31 + t27 * t51) + t53 * (t12 * t27 + t31 * t51) + t54 * t13) * g(1);
U = t1;
