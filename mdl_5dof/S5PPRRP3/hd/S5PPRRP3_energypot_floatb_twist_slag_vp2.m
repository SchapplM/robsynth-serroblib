% Calculate potential energy for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:29
% EndTime: 2019-12-05 15:10:29
% DurationCPUTime: 0.47s
% Computational Cost: add. (163->76), mult. (291->81), div. (0->0), fcn. (301->8), ass. (0->36)
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t58 = -t29 * mrSges(3,1) + t27 * mrSges(3,2) - mrSges(2,1);
t57 = mrSges(2,2) - mrSges(3,3);
t56 = -mrSges(5,3) - mrSges(6,2);
t55 = mrSges(4,2) + t56;
t54 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t53 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t52 = pkin(2) * t29;
t32 = sin(qJ(3));
t51 = t27 * t32;
t34 = cos(qJ(3));
t50 = t27 * t34;
t28 = sin(pkin(7));
t49 = t28 * t27;
t48 = t28 * t32;
t47 = t28 * t34;
t30 = cos(pkin(7));
t46 = t30 * t27;
t45 = t30 * t32;
t44 = t30 * t34;
t26 = qJ(1) + r_base(3);
t43 = t30 * pkin(1) + t28 * qJ(2) + r_base(1);
t41 = t28 * pkin(1) - t30 * qJ(2) + r_base(2);
t40 = pkin(5) * t46 + t30 * t52 + t43;
t39 = t27 * pkin(2) - t29 * pkin(5) + t26;
t38 = pkin(5) * t49 + t28 * t52 + t41;
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t10 = t29 * t44 + t48;
t9 = t29 * t45 - t47;
t8 = t29 * t47 - t45;
t7 = t29 * t48 + t44;
t1 = (-m(1) * r_base(3) - m(4) * t39 - mrSges(1,3) - mrSges(2,3) + t56 * t51 + t59 * (pkin(3) * t50 + pkin(6) * t51 + t39) + (-mrSges(3,2) + mrSges(4,3)) * t29 + (-t34 * mrSges(4,1) + t32 * mrSges(4,2) - mrSges(3,1)) * t27 + (-m(2) - m(3)) * t26 + t54 * (-t29 * t31 + t33 * t50) + t53 * (t29 * t33 + t31 * t50)) * g(3) + (-m(3) * t41 - m(4) * t38 - t8 * mrSges(4,1) - mrSges(4,3) * t49 - mrSges(1,2) + t60 * r_base(2) + t59 * (t8 * pkin(3) + t7 * pkin(6) + t38) - t57 * t30 + t58 * t28 + t54 * (t31 * t49 + t8 * t33) + t53 * (t8 * t31 - t33 * t49) + t55 * t7) * g(2) + (-m(3) * t43 - m(4) * t40 - t10 * mrSges(4,1) - mrSges(4,3) * t46 - mrSges(1,1) + t60 * r_base(1) + t59 * (t10 * pkin(3) + t9 * pkin(6) + t40) + t54 * (t10 * t33 + t31 * t46) + t58 * t30 + t53 * (t10 * t31 - t33 * t46) + t57 * t28 + t55 * t9) * g(1);
U = t1;
