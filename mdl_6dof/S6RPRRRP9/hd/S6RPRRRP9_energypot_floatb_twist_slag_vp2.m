% Calculate potential energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:03
% EndTime: 2019-03-09 06:26:03
% DurationCPUTime: 0.49s
% Computational Cost: add. (185->69), mult. (214->58), div. (0->0), fcn. (184->8), ass. (0->26)
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t7 = t22 * pkin(4) + pkin(3);
t18 = qJ(4) + qJ(5);
t9 = cos(t18);
t49 = -m(6) * t7 - m(7) * (pkin(5) * t9 + t7) - mrSges(4,1) - m(5) * pkin(3) - t22 * mrSges(5,1) + t19 * mrSges(5,2);
t25 = -pkin(9) - pkin(8);
t48 = m(6) * t25 + m(7) * (-qJ(6) + t25) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t47 = -m(1) - m(2);
t46 = -mrSges(6,1) - mrSges(7,1);
t45 = mrSges(6,2) + mrSges(7,2);
t44 = -m(4) - m(5) - m(6) - m(7);
t40 = pkin(4) * t19;
t8 = sin(t18);
t42 = -m(6) * t40 - m(7) * (pkin(5) * t8 + t40) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,1) * t19 - t22 * mrSges(5,2);
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t41 = t49 * t20 - t48 * t23 + mrSges(2,2) - mrSges(3,3);
t21 = sin(qJ(1));
t39 = t20 * t21;
t24 = cos(qJ(1));
t38 = t20 * t24;
t17 = pkin(6) + r_base(3);
t35 = t21 * pkin(1) + r_base(2);
t32 = t24 * pkin(1) + t21 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t17 + t44 * (pkin(2) + t17) + (t45 * t8 + t46 * t9 + t49) * t23 + t48 * t20) * g(3) + (-m(3) * t35 - mrSges(1,2) + t47 * r_base(2) + t46 * (t21 * t8 - t9 * t38) - t45 * (t21 * t9 + t8 * t38) + t44 * (t21 * pkin(7) + t35) + t42 * t21 + ((m(3) - t44) * qJ(2) - t41) * t24) * g(2) + (-m(3) * t32 - mrSges(1,1) + t47 * r_base(1) + t46 * (t24 * t8 + t9 * t39) - t45 * (t24 * t9 - t8 * t39) + t44 * (t24 * pkin(7) + t32) + t42 * t24 + t41 * t21) * g(1);
U  = t1;
