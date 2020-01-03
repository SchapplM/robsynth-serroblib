% Calculate potential energy for
% S5RRPRR10
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:41
% EndTime: 2019-12-31 20:23:42
% DurationCPUTime: 0.60s
% Computational Cost: add. (271->86), mult. (561->105), div. (0->0), fcn. (668->12), ass. (0->46)
t68 = -m(5) - m(6);
t67 = -mrSges(3,3) - mrSges(4,3);
t66 = -m(1) - m(2) - m(3);
t65 = -m(3) * pkin(1) - mrSges(2,1);
t64 = m(3) * pkin(7) - t67;
t63 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t62 = t35 * mrSges(6,1) + t39 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t61 = -m(6) * pkin(4) - mrSges(6,1) * t39 + mrSges(6,2) * t35 - mrSges(5,1);
t32 = sin(pkin(5));
t37 = sin(qJ(2));
t60 = t32 * t37;
t38 = sin(qJ(1));
t59 = t32 * t38;
t42 = cos(qJ(1));
t58 = t32 * t42;
t33 = cos(pkin(10));
t41 = cos(qJ(2));
t57 = t33 * t41;
t56 = t37 * t42;
t55 = t38 * t37;
t54 = t38 * t41;
t53 = t41 * t42;
t52 = pkin(6) + r_base(3);
t34 = cos(pkin(5));
t51 = t34 * pkin(7) + t52;
t19 = pkin(2) * t34 * t37 + (-pkin(7) - qJ(3)) * t32;
t28 = pkin(2) * t41 + pkin(1);
t50 = t42 * t19 + t38 * t28 + r_base(2);
t31 = sin(pkin(10));
t49 = t31 * t41 + t33 * t37;
t21 = -t31 * t37 + t57;
t48 = -t19 * t38 + t42 * t28 + r_base(1);
t47 = pkin(2) * t60 + t34 * qJ(3) + t51;
t46 = t21 * t34;
t40 = cos(qJ(4));
t36 = sin(qJ(4));
t18 = t49 * t34;
t17 = t49 * t32;
t16 = t31 * t60 - t32 * t57;
t10 = -t38 * t18 + t21 * t42;
t9 = -t38 * t46 - t42 * t49;
t8 = t18 * t42 + t38 * t21;
t7 = -t38 * t49 + t42 * t46;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t52 - mrSges(2,3) - m(3) * t51 - (mrSges(3,1) * t37 + mrSges(3,2) * t41) * t32 - m(4) * t47 - t17 * mrSges(4,1) + t68 * (t17 * pkin(3) + pkin(8) * t16 + t47) + t67 * t34 + t61 * (t17 * t40 + t34 * t36) - t62 * t16 + t63 * (t17 * t36 - t34 * t40)) * g(3) + (-mrSges(1,2) - t42 * mrSges(2,2) - (t34 * t56 + t54) * mrSges(3,1) - (t34 * t53 - t55) * mrSges(3,2) - m(4) * t50 - t8 * mrSges(4,1) + t65 * t38 + t66 * r_base(2) + t68 * (t8 * pkin(3) - pkin(8) * t7 + t50) + t61 * (-t36 * t58 + t8 * t40) + t62 * t7 + t64 * t58 + t63 * (t8 * t36 + t40 * t58)) * g(2) + (-mrSges(1,1) + t38 * mrSges(2,2) - (-t34 * t55 + t53) * mrSges(3,1) - (-t34 * t54 - t56) * mrSges(3,2) - m(4) * t48 - t10 * mrSges(4,1) + t65 * t42 + t66 * r_base(1) + t68 * (t10 * pkin(3) - pkin(8) * t9 + t48) + t61 * (t10 * t40 + t36 * t59) + t62 * t9 - t64 * t59 + t63 * (t10 * t36 - t40 * t59)) * g(1);
U = t1;
