% Calculate potential energy for
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:28
% EndTime: 2019-03-09 03:27:29
% DurationCPUTime: 0.58s
% Computational Cost: add. (192->74), mult. (229->61), div. (0->0), fcn. (203->8), ass. (0->33)
t21 = sin(pkin(9));
t22 = cos(pkin(9));
t59 = -m(5) * pkin(3) - t22 * mrSges(5,1) + t21 * mrSges(5,2) - mrSges(4,1);
t58 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t11 = pkin(4) * t22 + pkin(3);
t23 = -pkin(8) - qJ(4);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t57 = t11 * t24 + t23 * t26;
t56 = -m(1) - m(2);
t55 = -m(4) - m(5);
t54 = -m(6) - m(7);
t52 = -t21 * mrSges(5,1) - t22 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t51 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t50 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t49 = t59 * t24 - t58 * t26 + mrSges(2,2) - mrSges(3,3);
t48 = pkin(4) * t21;
t19 = pkin(9) + qJ(5);
t12 = sin(t19);
t27 = cos(qJ(1));
t46 = t12 * t27;
t25 = sin(qJ(1));
t44 = t25 * t12;
t13 = cos(t19);
t43 = t25 * t13;
t40 = t27 * t13;
t20 = pkin(6) + r_base(3);
t39 = t25 * pkin(1) + r_base(2);
t38 = pkin(2) + t20;
t37 = t25 * pkin(7) + t39;
t36 = t27 * pkin(1) + t25 * qJ(2) + r_base(1);
t34 = t27 * pkin(7) + t36;
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t55 * t38 + t54 * (t26 * t11 - t23 * t24 + t38) + (-m(2) - m(3)) * t20 + (-t50 * t12 + t51 * t13 + t59) * t26 + t58 * t24) * g(3) + (-m(3) * t39 - mrSges(1,2) + t56 * r_base(2) + t54 * (t25 * t48 + t37) + t51 * (-t24 * t40 + t44) + t55 * t37 + t50 * (t24 * t46 + t43) + t52 * t25 + (t54 * (-qJ(2) - t57) + (m(3) - t55) * qJ(2) - t49) * t27) * g(2) + (-m(3) * t36 - mrSges(1,1) + t56 * r_base(1) + t55 * t34 + t54 * (t27 * t48 + t34) + t51 * (t24 * t43 + t46) - t50 * (t24 * t44 - t40) + t52 * t27 + (t54 * t57 + t49) * t25) * g(1);
U  = t1;
