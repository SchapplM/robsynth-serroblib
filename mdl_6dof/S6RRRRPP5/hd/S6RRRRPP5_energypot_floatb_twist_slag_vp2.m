% Calculate potential energy for
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:00
% EndTime: 2019-03-09 21:05:01
% DurationCPUTime: 0.64s
% Computational Cost: add. (245->86), mult. (294->74), div. (0->0), fcn. (278->8), ass. (0->39)
t66 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3) + mrSges(7,3);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t65 = -m(4) * pkin(2) - mrSges(4,1) * t26 + mrSges(4,2) * t23 - mrSges(3,1);
t63 = -m(6) - m(7);
t24 = sin(qJ(2));
t27 = cos(qJ(2));
t29 = -pkin(9) - pkin(8);
t62 = t65 * t27 - mrSges(2,1) + (-m(7) * (-qJ(6) - t29) + t66) * t24;
t61 = -m(1) - m(2);
t60 = -m(3) - m(4);
t22 = qJ(3) + qJ(4);
t16 = sin(t22);
t17 = cos(t22);
t28 = cos(qJ(1));
t25 = sin(qJ(1));
t46 = t25 * t27;
t3 = t16 * t46 + t17 * t28;
t4 = -t16 * t28 + t17 * t46;
t59 = t4 * pkin(4) + t3 * qJ(5);
t45 = t27 * t28;
t5 = t16 * t45 - t25 * t17;
t6 = t16 * t25 + t17 * t45;
t58 = t6 * pkin(4) + t5 * qJ(5);
t56 = -t23 * mrSges(4,1) - t26 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t55 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t51 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t50 = pkin(3) * t23;
t47 = t24 * t29;
t21 = pkin(6) + r_base(3);
t42 = t25 * pkin(1) + r_base(2);
t41 = t28 * pkin(1) + t25 * pkin(7) + r_base(1);
t37 = -t28 * pkin(7) + t42;
t14 = pkin(3) * t26 + pkin(2);
t36 = t14 * t45 + t25 * t50 + t41;
t32 = -t28 * t47 + t36;
t7 = t14 * t46;
t31 = -t25 * t47 + t7 + (-pkin(7) - t50) * t28 + t42;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t60) * t21 + (-m(7) * qJ(6) - t66) * t27 + (-m(5) + t63) * (t24 * t14 + t27 * t29 + t21) + (t63 * (pkin(4) * t17 + qJ(5) * t16) + t55 * t16 + t51 * t17 + t65) * t24) * g(3) + (-mrSges(1,2) - m(3) * t37 - m(4) * t42 - m(5) * t31 - m(6) * (t31 + t59) - m(7) * (t37 + t7 + t59) + t61 * r_base(2) + t51 * t4 + t55 * t3 + (m(4) * pkin(7) + m(7) * t50 - t56) * t28 + t62 * t25) * g(2) + (-mrSges(1,1) - m(5) * t32 - m(6) * (t32 + t58) - m(7) * (t36 + t58) + t61 * r_base(1) + t60 * t41 + t51 * t6 + t55 * t5 + t56 * t25 + t62 * t28) * g(1);
U  = t1;
