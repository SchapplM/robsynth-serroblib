% Calculate potential energy for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:57
% EndTime: 2019-12-31 17:37:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (138->59), mult. (227->56), div. (0->0), fcn. (216->8), ass. (0->28)
t57 = -mrSges(3,1) - mrSges(4,1);
t56 = mrSges(3,2) - mrSges(4,3);
t27 = sin(qJ(2));
t29 = cos(qJ(2));
t26 = sin(qJ(5));
t28 = cos(qJ(5));
t47 = -m(6) * pkin(4) - t28 * mrSges(6,1) + t26 * mrSges(6,2) - mrSges(5,1);
t49 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t22 = sin(pkin(8));
t24 = cos(pkin(8));
t5 = t27 * t22 + t29 * t24;
t51 = t29 * t22 - t27 * t24;
t55 = t56 * t27 + t57 * t29 + t5 * t47 + t51 * t49 - mrSges(2,1);
t52 = -m(5) - m(6);
t54 = -t26 * mrSges(6,1) - t28 * mrSges(6,2) + t52 * qJ(4) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t53 = -m(1) - m(2);
t23 = sin(pkin(7));
t44 = t23 * t29;
t25 = cos(pkin(7));
t43 = t25 * t29;
t40 = qJ(3) * t27;
t21 = qJ(1) + r_base(3);
t39 = t25 * pkin(1) + t23 * pkin(5) + r_base(1);
t36 = t23 * pkin(1) - t25 * pkin(5) + r_base(2);
t35 = pkin(2) * t43 + t25 * t40 + t39;
t34 = t27 * pkin(2) - t29 * qJ(3) + t21;
t33 = pkin(2) * t44 + t23 * t40 + t36;
t1 = (-m(1) * r_base(3) - m(4) * t34 - mrSges(1,3) - mrSges(2,3) - t47 * t51 + t52 * (t27 * pkin(3) + t34) + t49 * t5 - t56 * t29 + t57 * t27 + (-m(2) - m(3)) * t21) * g(3) + (-m(3) * t36 - m(4) * t33 - mrSges(1,2) + t53 * r_base(2) + t52 * (pkin(3) * t44 + t33)) * g(2) + (-m(3) * t39 - m(4) * t35 - mrSges(1,1) + t53 * r_base(1) + t52 * (pkin(3) * t43 + t35)) * g(1) + (t55 * g(1) + t54 * g(2)) * t25 + (-t54 * g(1) + t55 * g(2)) * t23;
U = t1;
