% Calculate potential energy for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR14_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:41
% EndTime: 2019-12-31 20:35:41
% DurationCPUTime: 0.71s
% Computational Cost: add. (245->98), mult. (383->110), div. (0->0), fcn. (422->12), ass. (0->43)
t64 = -m(1) - m(2);
t63 = -m(5) - m(6);
t62 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t61 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t60 = -m(6) * pkin(4) - t38 * mrSges(6,1) + t35 * mrSges(6,2) - mrSges(5,1);
t59 = -t35 * mrSges(6,1) - t38 * mrSges(6,2) - mrSges(5,3) - t61;
t30 = sin(pkin(10));
t58 = pkin(3) * t30;
t31 = sin(pkin(5));
t36 = sin(qJ(2));
t57 = t31 * t36;
t37 = sin(qJ(1));
t56 = t31 * t37;
t39 = cos(qJ(2));
t55 = t31 * t39;
t40 = cos(qJ(1));
t54 = t31 * t40;
t53 = t36 * t37;
t52 = t36 * t40;
t51 = t37 * t39;
t50 = t39 * t40;
t49 = pkin(6) + r_base(3);
t48 = t37 * pkin(1) + r_base(2);
t47 = t30 * t56;
t33 = cos(pkin(5));
t46 = t33 * pkin(7) + t49;
t45 = t40 * pkin(1) + pkin(7) * t56 + r_base(1);
t44 = -pkin(7) * t54 + t48;
t32 = cos(pkin(10));
t23 = pkin(3) * t32 + pkin(2);
t34 = -pkin(8) - qJ(3);
t43 = t23 * t57 + t33 * t58 + t34 * t55 + t46;
t29 = pkin(10) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t14 = -t33 * t53 + t50;
t13 = t33 * t51 + t52;
t12 = t33 * t52 + t51;
t11 = -t33 * t50 + t53;
t8 = t24 * t33 + t25 * t57;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t49 - mrSges(2,3) - m(5) * t43 - t8 * mrSges(5,1) + mrSges(5,3) * t55 - m(6) * (pkin(4) * t8 + t43) - (-t35 * t55 + t38 * t8) * mrSges(6,1) - (-t35 * t8 - t38 * t55) * mrSges(6,2) + t62 * (t24 * t57 - t33 * t25) + (-m(3) - m(4)) * t46 + (-t30 * mrSges(4,1) - t32 * mrSges(4,2) - mrSges(3,3)) * t33 + (t61 * t39 + (-m(4) * pkin(2) - t32 * mrSges(4,1) + t30 * mrSges(4,2) - mrSges(3,1)) * t36) * t31) * g(3) + (-mrSges(1,2) - t37 * mrSges(2,1) - t40 * mrSges(2,2) - m(3) * t44 - t12 * mrSges(3,1) + mrSges(3,3) * t54 - m(4) * (t12 * pkin(2) + t44) - (t12 * t32 - t30 * t54) * mrSges(4,1) - (-t12 * t30 - t32 * t54) * mrSges(4,2) + t64 * r_base(2) + t63 * (t12 * t23 - t11 * t34 + (-pkin(7) - t58) * t54 + t48) + t62 * (t12 * t24 + t25 * t54) + t60 * (t12 * t25 - t24 * t54) + t59 * t11) * g(2) + (-mrSges(1,1) - t40 * mrSges(2,1) + t37 * mrSges(2,2) - m(3) * t45 - t14 * mrSges(3,1) - mrSges(3,3) * t56 - m(4) * (pkin(2) * t14 + t45) - (t14 * t32 + t47) * mrSges(4,1) - (-t14 * t30 + t32 * t56) * mrSges(4,2) + t64 * r_base(1) + t63 * (pkin(3) * t47 - t13 * t34 + t14 * t23 + t45) + t62 * (t14 * t24 - t25 * t56) + t60 * (t14 * t25 + t24 * t56) + t59 * t13) * g(1);
U = t1;
