% Calculate potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:52
% DurationCPUTime: 0.55s
% Computational Cost: add. (151->78), mult. (216->68), div. (0->0), fcn. (186->6), ass. (0->28)
t51 = -mrSges(4,1) + mrSges(5,2);
t50 = -m(1) - m(2);
t49 = -m(6) - m(7);
t48 = -mrSges(6,3) - mrSges(7,2);
t47 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t46 = -t24 * mrSges(4,2) + t51 * t21 + mrSges(2,2) - mrSges(3,3);
t45 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t44 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t22 = sin(qJ(1));
t41 = t21 * t22;
t40 = t22 * t24;
t25 = cos(qJ(1));
t39 = t24 * t25;
t38 = qJ(4) * t24;
t19 = pkin(6) + r_base(3);
t37 = t22 * pkin(1) + r_base(2);
t36 = pkin(2) + t19;
t35 = t22 * pkin(7) + t37;
t34 = t25 * pkin(1) + t22 * qJ(2) + r_base(1);
t33 = t25 * t38 + t35;
t32 = t25 * pkin(7) + t34;
t31 = t24 * pkin(3) + t21 * qJ(4) + t36;
t29 = pkin(3) * t41 + t32;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - m(4) * t36 - m(5) * t31 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t49 * (t24 * pkin(8) + t31) + (-m(2) - m(3)) * t19 + (t48 + t51) * t24 + (t45 * t20 + t44 * t23 + mrSges(4,2) - mrSges(5,3)) * t21) * g(3) + (-m(3) * t37 - m(4) * t35 - m(5) * t33 - mrSges(1,2) + t50 * r_base(2) + t49 * (t22 * pkin(4) + t33) + t45 * (t20 * t39 + t22 * t23) - t44 * (t20 * t22 - t23 * t39) + t47 * t22 + (-t24 * mrSges(5,3) + (m(5) * pkin(3) + t49 * (-pkin(3) - pkin(8)) - t48) * t21 + (m(3) + m(4) + m(5) - t49) * qJ(2) - t46) * t25) * g(2) + (-m(3) * t34 - m(4) * t32 - m(5) * t29 - mrSges(1,1) + t50 * r_base(1) + t48 * t41 + t49 * (t25 * pkin(4) + pkin(8) * t41 - t22 * t38 + t29) + t45 * (-t20 * t40 + t23 * t25) - t44 * (t20 * t25 + t23 * t40) + t47 * t25 + (-(-m(5) * qJ(4) - mrSges(5,3)) * t24 + t46) * t22) * g(1);
U  = t1;
