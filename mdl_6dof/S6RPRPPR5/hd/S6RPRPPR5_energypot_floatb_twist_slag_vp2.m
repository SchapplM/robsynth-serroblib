% Calculate potential energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:02
% EndTime: 2019-03-09 02:50:02
% DurationCPUTime: 0.56s
% Computational Cost: add. (227->83), mult. (218->72), div. (0->0), fcn. (184->10), ass. (0->37)
t16 = pkin(10) + qJ(6);
t11 = sin(t16);
t13 = cos(t16);
t19 = sin(pkin(10));
t60 = -m(7) * pkin(5) * t19 - t11 * mrSges(7,1) - t13 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t59 = -mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5);
t57 = -m(2) - m(3);
t56 = -m(6) - m(7);
t26 = cos(qJ(1));
t17 = pkin(9) + qJ(3);
t12 = sin(t17);
t41 = qJ(4) * t12;
t14 = cos(t17);
t46 = t14 * t26;
t55 = pkin(3) * t46 + t26 * t41;
t54 = -m(1) + t57;
t53 = m(5) - t56;
t20 = sin(pkin(9));
t22 = cos(pkin(9));
t50 = -m(3) * pkin(1) - t22 * mrSges(3,1) + t20 * mrSges(3,2) + t60 * t12 + t59 * t14 - mrSges(2,1);
t49 = -m(3) * qJ(2) - t13 * mrSges(7,1) + t11 * mrSges(7,2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t25 = sin(qJ(1));
t47 = t14 * t25;
t45 = t19 * t26;
t21 = cos(pkin(10));
t44 = t21 * t26;
t43 = t25 * t19;
t42 = t25 * t21;
t18 = pkin(6) + r_base(3);
t9 = pkin(2) * t22 + pkin(1);
t39 = t26 * t9 + r_base(1);
t38 = t20 * pkin(2) + t18;
t24 = -pkin(7) - qJ(2);
t37 = t26 * t24 + t25 * t9 + r_base(2);
t29 = -t25 * t24 + t39;
t8 = pkin(5) * t21 + pkin(4);
t1 = (-m(1) * r_base(3) - m(4) * t38 - t20 * mrSges(3,1) - t22 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t57 * t18 - t53 * (t12 * pkin(3) + t38) + (t19 * mrSges(6,1) + t21 * mrSges(6,2) + qJ(4) * t53 - t60) * t14 + (-mrSges(6,3) + t59) * t12) * g(3) + (-mrSges(1,2) - m(4) * t37 - (t12 * t43 - t44) * mrSges(6,1) - (t12 * t42 + t45) * mrSges(6,2) - mrSges(6,3) * t47 + t54 * r_base(2) - t53 * (pkin(3) * t47 + t25 * t41 + t37) + (m(6) * pkin(4) + m(7) * t8 - t49) * t26 + t50 * t25) * g(2) + (-mrSges(1,1) - m(4) * t29 - m(5) * (t29 + t55) - (t12 * t45 + t42) * mrSges(6,1) - (t12 * t44 - t43) * mrSges(6,2) - mrSges(6,3) * t46 + t56 * (t39 + t55) + t54 * r_base(1) + t50 * t26 + (-m(6) * (pkin(4) - t24) - m(7) * (-t24 + t8) + t49) * t25) * g(1);
U  = t1;
