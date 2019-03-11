% Calculate potential energy for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:09
% EndTime: 2019-03-09 09:34:10
% DurationCPUTime: 0.57s
% Computational Cost: add. (207->85), mult. (251->73), div. (0->0), fcn. (221->10), ass. (0->35)
t19 = pkin(10) + qJ(5);
t10 = sin(t19);
t11 = cos(t19);
t21 = sin(pkin(10));
t49 = t21 * pkin(4);
t12 = qJ(6) + t19;
t7 = sin(t12);
t8 = cos(t12);
t62 = -m(6) * t49 - t10 * mrSges(6,1) - t11 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t10 + t49) - t7 * mrSges(7,1) - t8 * mrSges(7,2);
t23 = -pkin(8) - qJ(4);
t61 = -mrSges(3,1) + mrSges(4,2) + m(7) * (-pkin(9) + t23) - mrSges(7,3) + m(6) * t23 - mrSges(6,3) - m(5) * qJ(4);
t59 = -m(1) - m(2);
t25 = sin(qJ(1));
t24 = sin(qJ(2));
t43 = qJ(3) * t24;
t26 = cos(qJ(2));
t45 = t25 * t26;
t58 = pkin(2) * t45 + t25 * t43;
t57 = -m(5) - m(6) - m(7);
t56 = m(4) - t57;
t51 = t11 * mrSges(6,1) + t8 * mrSges(7,1) - t10 * mrSges(6,2) - t7 * mrSges(7,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t50 = t62 * t24 + t61 * t26 - mrSges(2,1);
t22 = cos(pkin(10));
t9 = t22 * pkin(4) + pkin(3);
t27 = cos(qJ(1));
t48 = t24 * t27;
t47 = t25 * t21;
t46 = t25 * t22;
t44 = t26 * t27;
t20 = pkin(6) + r_base(3);
t41 = t25 * pkin(1) + r_base(2);
t39 = t27 * pkin(1) + t25 * pkin(7) + r_base(1);
t31 = -pkin(7) * t27 + t41;
t1 = pkin(5) * t11 + t9;
t2 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t20 - t56 * (t24 * pkin(2) + t20) + (t21 * mrSges(5,1) + t22 * mrSges(5,2) + t56 * qJ(3) - t62) * t26 + (-mrSges(5,3) + t61) * t24) * g(3) + (-mrSges(1,2) - m(3) * t31 - m(4) * (t31 + t58) - mrSges(5,3) * t45 + t59 * r_base(2) + (-mrSges(5,1) * t47 - mrSges(5,2) * t46) * t24 + t57 * (t41 + t58) + (-m(5) * (-pkin(3) - pkin(7)) + t22 * mrSges(5,1) - t21 * mrSges(5,2) - m(6) * (-pkin(7) - t9) - m(7) * (-pkin(7) - t1) + t51) * t27 + t50 * t25) * g(2) + (-mrSges(1,1) - m(3) * t39 - (t21 * t48 + t46) * mrSges(5,1) - (t22 * t48 - t47) * mrSges(5,2) - mrSges(5,3) * t44 + t59 * r_base(1) - t56 * (pkin(2) * t44 + t27 * t43 + t39) + t50 * t27 + (-m(5) * pkin(3) - m(6) * t9 - m(7) * t1 - t51) * t25) * g(1);
U  = t2;
