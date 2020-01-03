% Calculate potential energy for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:06
% EndTime: 2019-12-31 17:18:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (92->48), mult. (135->45), div. (0->0), fcn. (113->6), ass. (0->20)
t14 = cos(qJ(3));
t34 = -m(4) * pkin(2) - m(5) * (pkin(3) * t14 + pkin(2)) - mrSges(3,1);
t33 = m(4) * pkin(6) - m(5) * (-qJ(4) - pkin(6)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t32 = m(5) * pkin(3);
t31 = -m(1) - m(2);
t30 = -mrSges(4,1) - mrSges(5,1);
t29 = mrSges(2,2) - mrSges(3,3);
t28 = mrSges(4,2) + mrSges(5,2);
t27 = -m(3) - m(4) - m(5);
t12 = sin(qJ(2));
t15 = cos(qJ(2));
t26 = -t33 * t12 + t34 * t15 - mrSges(2,1);
t11 = sin(qJ(3));
t16 = cos(qJ(1));
t25 = t11 * t16;
t13 = sin(qJ(1));
t24 = t13 * t11;
t23 = t13 * t15;
t22 = t16 * t15;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t27) * (pkin(4) + r_base(3)) + t33 * t15 + (t28 * t11 + t30 * t14 + t34) * t12) * g(3) + (t25 * t32 - mrSges(1,2) + t31 * r_base(2) + t27 * (t13 * pkin(1) - pkin(5) * t16 + r_base(2)) + t30 * (t14 * t23 - t25) - t29 * t16 - t28 * (-t11 * t23 - t14 * t16) + t26 * t13) * g(2) + (-t24 * t32 - mrSges(1,1) + t31 * r_base(1) + t30 * (t14 * t22 + t24) - t28 * (-t11 * t22 + t13 * t14) + t27 * (t16 * pkin(1) + t13 * pkin(5) + r_base(1)) + t29 * t13 + t26 * t16) * g(1);
U = t1;
