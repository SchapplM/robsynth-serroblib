% Calculate potential energy for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:56
% EndTime: 2019-12-31 20:28:56
% DurationCPUTime: 0.47s
% Computational Cost: add. (138->62), mult. (227->59), div. (0->0), fcn. (216->8), ass. (0->29)
t55 = -mrSges(3,1) - mrSges(4,1);
t54 = mrSges(3,2) - mrSges(4,3);
t24 = sin(qJ(2));
t28 = cos(qJ(2));
t22 = sin(qJ(5));
t26 = cos(qJ(5));
t47 = -m(6) * pkin(4) - t26 * mrSges(6,1) + t22 * mrSges(6,2) - mrSges(5,1);
t23 = sin(qJ(4));
t27 = cos(qJ(4));
t5 = t24 * t23 + t28 * t27;
t53 = t54 * t24 + t55 * t28 + t5 * t47 - mrSges(2,1);
t52 = -m(1) - m(2);
t51 = -m(5) - m(6);
t44 = t24 * t27;
t50 = t28 * t23 - t44;
t48 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t46 = t22 * mrSges(6,1) + t26 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t25 = sin(qJ(1));
t43 = t25 * t28;
t29 = cos(qJ(1));
t41 = t29 * t28;
t40 = qJ(3) * t24;
t21 = pkin(5) + r_base(3);
t39 = t29 * pkin(1) + t25 * pkin(6) + r_base(1);
t36 = t25 * pkin(1) - t29 * pkin(6) + r_base(2);
t35 = pkin(2) * t41 + t29 * t40 + t39;
t34 = t24 * pkin(2) - t28 * qJ(3) + t21;
t33 = pkin(2) * t43 + t25 * t40 + t36;
t1 = (-m(1) * r_base(3) - m(4) * t34 - mrSges(1,3) - mrSges(2,3) - t47 * t50 + t51 * (t24 * pkin(3) + t34) + t48 * t5 - t54 * t28 + t55 * t24 + (-m(2) - m(3)) * t21) * g(3) + (-m(3) * t36 - m(4) * t33 - mrSges(1,2) + t52 * r_base(2) + t51 * (pkin(3) * t43 + t29 * pkin(7) + t33) - t46 * t29) * g(2) + (-m(3) * t39 - m(4) * t35 - mrSges(1,1) + t52 * r_base(1) + t51 * (pkin(3) * t41 + t35) + t48 * t23 * t41 + (-t48 * t44 + t53) * t29) * g(1) + ((t48 * t50 + t53) * g(2) + (-t51 * pkin(7) + t46) * g(1)) * t25;
U = t1;
