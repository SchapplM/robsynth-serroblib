% Calculate potential energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:13
% EndTime: 2019-03-08 18:12:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (75->41), mult. (90->32), div. (0->0), fcn. (72->4), ass. (0->17)
t27 = -m(4) - m(5);
t13 = sin(pkin(5));
t14 = cos(pkin(5));
t25 = t14 * pkin(1) + t13 * qJ(2);
t24 = cos(qJ(3));
t23 = sin(qJ(3));
t22 = -mrSges(2,1) - mrSges(3,1);
t21 = mrSges(2,2) - mrSges(3,3);
t20 = m(3) - t27;
t18 = qJ(1) + r_base(3);
t17 = -m(1) - m(2) - t20;
t16 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t15 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t9 = t13 * pkin(1);
t2 = -t13 * t24 + t14 * t23;
t1 = -t13 * t23 - t14 * t24;
t3 = (-m(1) * r_base(3) - mrSges(3,2) + mrSges(5,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t27 * (-pkin(4) + t18) + (-m(2) - m(3)) * t18) * g(3) + (-mrSges(1,2) - m(3) * t9 + t17 * r_base(2) + t16 * t2 + (t20 * qJ(2) - t21) * t14 + t22 * t13 + t15 * t1 + t27 * (t13 * pkin(2) + t9)) * g(2) + (-m(3) * t25 + t16 * t1 + t21 * t13 + t22 * t14 - t15 * t2 + t17 * r_base(1) - mrSges(1,1) + t27 * (t14 * pkin(2) + t25)) * g(1);
U  = t3;
