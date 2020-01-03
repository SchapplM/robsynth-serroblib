% Calculate potential energy for
% S4RRPP5
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:07
% EndTime: 2019-12-31 17:00:07
% DurationCPUTime: 0.28s
% Computational Cost: add. (79->43), mult. (109->33), div. (0->0), fcn. (79->4), ass. (0->16)
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t32 = pkin(2) * t12 + qJ(3) * t10;
t31 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t30 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t29 = -m(1) - m(2);
t28 = -m(4) - m(5);
t11 = sin(qJ(1));
t27 = t32 * t11;
t25 = t30 * t10 + t31 * t12 - mrSges(2,1);
t24 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3) - mrSges(2,2);
t9 = pkin(4) + r_base(3);
t21 = t11 * pkin(1) + r_base(2);
t13 = cos(qJ(1));
t18 = -pkin(5) * t13 + t21;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 + t28 * (t10 * pkin(2) - qJ(3) * t12 + t9) - t30 * t12 + t31 * t10) * g(3) + (-mrSges(1,2) - m(3) * t18 - m(4) * (t18 + t27) - m(5) * (t21 + t27) + t29 * r_base(2) + (-m(5) * (-pkin(3) - pkin(5)) + t24) * t13 + t25 * t11) * g(2) + (-mrSges(1,1) + t29 * r_base(1) + (-m(5) * pkin(3) - t24) * t11 + (-m(3) + t28) * (t13 * pkin(1) + t11 * pkin(5) + r_base(1)) + (t28 * t32 + t25) * t13) * g(1);
U = t1;
