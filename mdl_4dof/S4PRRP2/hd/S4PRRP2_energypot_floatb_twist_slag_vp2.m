% Calculate potential energy for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:58
% EndTime: 2019-03-08 18:23:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (64->38), mult. (48->24), div. (0->0), fcn. (18->4), ass. (0->13)
t17 = -m(4) - m(5);
t16 = -m(2) - m(3);
t15 = pkin(5) + pkin(4);
t13 = mrSges(4,2) + mrSges(5,2);
t12 = qJ(1) + r_base(2);
t10 = -m(1) + t16 + t17;
t9 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t8 = cos(qJ(2));
t7 = sin(qJ(2));
t6 = qJ(2) + qJ(3);
t2 = cos(t6);
t1 = sin(t6);
t3 = (-mrSges(1,3) + mrSges(2,2) - m(3) * pkin(4) - mrSges(3,3) - m(4) * t15 - mrSges(4,3) - m(5) * (qJ(4) + t15) - mrSges(5,3) + t10 * r_base(3)) * g(3) + (-m(1) * r_base(2) - t7 * mrSges(3,1) - t8 * mrSges(3,2) + t9 * t1 + t16 * t12 - t13 * t2 - mrSges(1,2) - mrSges(2,3) + t17 * (t7 * pkin(2) + t12)) * g(2) + (-m(3) * pkin(1) - t8 * mrSges(3,1) + t7 * mrSges(3,2) + t13 * t1 + t10 * r_base(1) + t9 * t2 - mrSges(1,1) - mrSges(2,1) + t17 * (t8 * pkin(2) + pkin(1))) * g(1);
U  = t3;
