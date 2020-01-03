% Calculate potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:16
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.50s
% Computational Cost: add. (186->77), mult. (379->83), div. (0->0), fcn. (406->8), ass. (0->44)
t68 = -m(1) - m(2);
t67 = -m(5) - m(6);
t29 = cos(pkin(8));
t27 = sin(pkin(8));
t31 = sin(qJ(2));
t54 = t31 * t27;
t30 = cos(pkin(5));
t33 = cos(qJ(2));
t55 = t30 * t33;
t8 = -t29 * t55 + t54;
t53 = t31 * t29;
t9 = t27 * t55 + t53;
t66 = t9 * pkin(3) + t8 * qJ(4);
t65 = -t33 * mrSges(3,1) + t31 * mrSges(3,2) - mrSges(2,1);
t64 = mrSges(2,2) - mrSges(3,3);
t32 = sin(qJ(1));
t52 = t32 * t31;
t28 = sin(pkin(5));
t34 = cos(qJ(1));
t56 = t28 * t34;
t63 = t30 * t52 + t56;
t62 = mrSges(5,1) + mrSges(6,1) + mrSges(4,3);
t61 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t60 = -m(6) * pkin(4) - t62;
t59 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t58 = t28 * t32;
t51 = t32 * t33;
t50 = t33 * t34;
t49 = t34 * t30;
t47 = qJ(3) * t28;
t46 = qJ(3) * t30;
t26 = pkin(6) + r_base(3);
t45 = t32 * pkin(1) + r_base(2);
t43 = t31 * t47;
t42 = t31 * pkin(2) + t26;
t41 = t34 * pkin(1) + t32 * pkin(7) + r_base(1);
t39 = pkin(2) * t50 + t32 * t46 + t34 * t43 + t41;
t38 = -t33 * t47 + t42;
t36 = t32 * t43 + pkin(2) * t51 + (-pkin(7) - t46) * t34 + t45;
t6 = t27 * t58 + t29 * t50 - t49 * t54;
t5 = -t29 * t58 + (t27 * t33 + t30 * t53) * t34;
t4 = -t27 * t63 + t29 * t51;
t3 = t27 * t51 + t29 * t63;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t31 * mrSges(3,1) - m(4) * t38 - m(5) * (t38 + t66) - m(6) * (t42 + t66) + (-m(2) - m(3)) * t26 + t59 * t9 + t61 * t8 + (-mrSges(3,2) + (-m(6) * (-pkin(4) - qJ(3)) + t62) * t28) * t33) * g(3) + (-m(3) * t45 - m(4) * t36 - mrSges(1,2) + t68 * r_base(2) + t67 * (t4 * pkin(3) + t3 * qJ(4) + t36) + (m(3) * pkin(7) - t64) * t34 + t65 * t32 + t59 * t4 + t61 * t3 + t60 * (t28 * t52 - t49)) * g(2) + (-m(3) * t41 - m(4) * t39 - mrSges(1,1) + t68 * r_base(1) + t67 * (t6 * pkin(3) + t5 * qJ(4) + t39) + t65 * t34 + t64 * t32 + t59 * t6 + t61 * t5 + t60 * (t32 * t30 + t31 * t56)) * g(1);
U = t1;
