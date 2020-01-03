% Calculate potential energy for
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:03
% EndTime: 2019-12-31 22:08:04
% DurationCPUTime: 0.54s
% Computational Cost: add. (225->92), mult. (447->106), div. (0->0), fcn. (511->10), ass. (0->44)
t63 = -m(1) - m(2);
t62 = -mrSges(5,1) - mrSges(6,1);
t61 = -mrSges(5,2) - mrSges(6,2);
t31 = sin(qJ(4));
t48 = pkin(4) * t31 + pkin(8);
t60 = -m(6) * t48 + mrSges(3,2) - mrSges(4,3);
t59 = -m(5) * pkin(9) + m(6) * (-qJ(5) - pkin(9)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t28 = sin(pkin(5));
t33 = sin(qJ(2));
t58 = t28 * t33;
t34 = sin(qJ(1));
t57 = t28 * t34;
t36 = cos(qJ(3));
t56 = t28 * t36;
t37 = cos(qJ(2));
t55 = t28 * t37;
t38 = cos(qJ(1));
t54 = t28 * t38;
t53 = t33 * t34;
t52 = t33 * t38;
t51 = t34 * t37;
t50 = t37 * t38;
t49 = pkin(6) + r_base(3);
t29 = cos(pkin(5));
t47 = t29 * pkin(7) + t49;
t46 = t38 * pkin(1) + pkin(7) * t57 + r_base(1);
t45 = pkin(2) * t58 + t47;
t18 = -t29 * t53 + t50;
t44 = t18 * pkin(2) + t46;
t43 = t34 * pkin(1) - pkin(7) * t54 + r_base(2);
t16 = t29 * t52 + t51;
t42 = t16 * pkin(2) + t43;
t17 = t29 * t51 + t52;
t41 = pkin(8) * t17 + t44;
t40 = -pkin(8) * t55 + t45;
t15 = -t29 * t50 + t53;
t39 = t15 * pkin(8) + t42;
t35 = cos(qJ(4));
t32 = sin(qJ(3));
t24 = pkin(4) * t35 + pkin(3);
t14 = t29 * t32 + t33 * t56;
t10 = t18 * t36 + t32 * t57;
t8 = t16 * t36 - t32 * t54;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t49 - mrSges(2,3) - m(3) * t47 - t29 * mrSges(3,3) - (t33 * mrSges(3,1) + t37 * mrSges(3,2)) * t28 - m(4) * t40 - t14 * mrSges(4,1) + mrSges(4,3) * t55 - m(5) * (pkin(3) * t14 + t40) - m(6) * (t14 * t24 - t48 * t55 + t45) + t62 * (t14 * t35 - t31 * t55) + t61 * (-t14 * t31 - t35 * t55) + t59 * (-t29 * t36 + t32 * t58)) * g(3) + (-mrSges(1,2) - t34 * mrSges(2,1) - t38 * mrSges(2,2) - m(3) * t43 - t16 * mrSges(3,1) + mrSges(3,3) * t54 - m(4) * t39 - t8 * mrSges(4,1) - m(5) * (t8 * pkin(3) + t39) - m(6) * (t8 * t24 + t42) + t63 * r_base(2) + t62 * (t15 * t31 + t35 * t8) + t60 * t15 + t61 * (t15 * t35 - t31 * t8) + t59 * (t16 * t32 + t36 * t54)) * g(2) + (-mrSges(1,1) - t38 * mrSges(2,1) + t34 * mrSges(2,2) - m(3) * t46 - t18 * mrSges(3,1) - mrSges(3,3) * t57 - m(4) * t41 - t10 * mrSges(4,1) - m(5) * (pkin(3) * t10 + t41) - m(6) * (t10 * t24 + t44) + t63 * r_base(1) + t62 * (t10 * t35 + t17 * t31) + t61 * (-t10 * t31 + t17 * t35) + t60 * t17 + t59 * (t18 * t32 - t34 * t56)) * g(1);
U = t1;
