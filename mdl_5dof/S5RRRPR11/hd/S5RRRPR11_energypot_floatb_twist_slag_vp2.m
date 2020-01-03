% Calculate potential energy for
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:20
% EndTime: 2019-12-31 21:32:20
% DurationCPUTime: 0.51s
% Computational Cost: add. (149->70), mult. (254->68), div. (0->0), fcn. (250->8), ass. (0->33)
t48 = -m(1) - m(2);
t47 = -m(5) - m(6);
t20 = sin(qJ(3));
t21 = sin(qJ(2));
t24 = cos(qJ(3));
t46 = (pkin(3) * t24 + qJ(4) * t20) * t21;
t25 = cos(qJ(2));
t45 = -t25 * mrSges(3,1) + t21 * mrSges(3,2) - mrSges(2,1);
t44 = mrSges(2,2) - mrSges(3,3);
t43 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t42 = m(6) * pkin(8) - t43;
t19 = sin(qJ(5));
t23 = cos(qJ(5));
t41 = -t19 * mrSges(6,1) - t23 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t40 = -m(6) * pkin(4) - t23 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(4,1) - mrSges(5,1);
t22 = sin(qJ(1));
t39 = t22 * t21;
t38 = t22 * t25;
t26 = cos(qJ(1));
t37 = t26 * t21;
t36 = t26 * t25;
t18 = pkin(5) + r_base(3);
t35 = t21 * pkin(2) + t18;
t34 = t26 * pkin(1) + t22 * pkin(6) + r_base(1);
t32 = t22 * pkin(1) - t26 * pkin(6) + r_base(2);
t31 = pkin(2) * t36 + pkin(7) * t37 + t34;
t30 = -t25 * pkin(7) + t35;
t29 = pkin(2) * t38 + pkin(7) * t39 + t32;
t6 = t22 * t20 + t24 * t36;
t5 = t20 * t36 - t22 * t24;
t4 = -t26 * t20 + t24 * t38;
t3 = t20 * t38 + t26 * t24;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t30 - m(5) * (t30 + t46) - m(6) * (t35 + t46) + (-m(2) - m(3)) * t18 + (-mrSges(3,2) - m(6) * (-pkin(7) + pkin(8)) + t43) * t25 + (t41 * t20 + t40 * t24 - mrSges(3,1)) * t21) * g(3) + (-m(3) * t32 - m(4) * t29 - mrSges(1,2) + t48 * r_base(2) + t47 * (t4 * pkin(3) + t3 * qJ(4) + t29) + t40 * t4 + t41 * t3 - t44 * t26 + t45 * t22 + t42 * t39) * g(2) + (-m(3) * t34 - m(4) * t31 - mrSges(1,1) + t48 * r_base(1) + t47 * (t6 * pkin(3) + t5 * qJ(4) + t31) + t40 * t6 + t41 * t5 + t45 * t26 + t44 * t22 + t42 * t37) * g(1);
U = t1;
