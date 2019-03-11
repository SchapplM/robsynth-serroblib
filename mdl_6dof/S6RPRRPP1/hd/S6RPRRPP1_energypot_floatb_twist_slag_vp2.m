% Calculate potential energy for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:48
% EndTime: 2019-03-09 04:27:49
% DurationCPUTime: 0.51s
% Computational Cost: add. (271->77), mult. (227->68), div. (0->0), fcn. (201->10), ass. (0->34)
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t58 = -m(5) * pkin(3) - t28 * mrSges(5,1) + t25 * mrSges(5,2) - mrSges(4,1);
t57 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t56 = -m(1) - m(2);
t55 = m(4) + m(5);
t54 = -m(6) - m(7);
t53 = -t25 * mrSges(5,1) - t28 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t51 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t50 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t49 = t57 * t26 + t58 * t29 - mrSges(3,1);
t48 = pkin(4) * t25;
t23 = qJ(1) + pkin(9);
t16 = sin(t23);
t47 = t16 * t29;
t18 = cos(t23);
t46 = t18 * t29;
t24 = -qJ(5) - pkin(8);
t45 = t24 * t26;
t42 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t41 = t27 * pkin(1) + r_base(2);
t30 = cos(qJ(1));
t40 = t30 * pkin(1) + r_base(1);
t19 = qJ(2) + t42;
t39 = t16 * pkin(2) + t41;
t38 = t18 * pkin(2) + t16 * pkin(7) + t40;
t22 = qJ(4) + pkin(10);
t17 = cos(t22);
t15 = sin(t22);
t14 = pkin(4) * t28 + pkin(3);
t1 = (-m(1) * r_base(3) - m(2) * t42 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t54 * (t26 * t14 + t29 * t24 + t19) + (-m(3) - t55) * t19 - t57 * t29 + (t15 * t50 + t17 * t51 + t58) * t26) * g(3) + (-m(3) * t41 - t27 * mrSges(2,1) - t30 * mrSges(2,2) - mrSges(1,2) + t56 * r_base(2) - t55 * t39 + t54 * (-t16 * t45 + t14 * t47 + (-pkin(7) - t48) * t18 + t39) + t51 * (-t15 * t18 + t17 * t47) + t50 * (t15 * t47 + t17 * t18) + (t55 * pkin(7) - t53) * t18 + t49 * t16) * g(2) + (-m(3) * t40 - t30 * mrSges(2,1) + t27 * mrSges(2,2) - mrSges(1,1) + t56 * r_base(1) + t54 * (t14 * t46 + t16 * t48 - t18 * t45 + t38) + t51 * (t15 * t16 + t17 * t46) - t55 * t38 + t50 * (t15 * t46 - t16 * t17) + t53 * t16 + t49 * t18) * g(1);
U  = t1;
