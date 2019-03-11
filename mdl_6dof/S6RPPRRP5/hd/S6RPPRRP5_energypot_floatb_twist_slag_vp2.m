% Calculate potential energy for
% S6RPPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:37
% EndTime: 2019-03-09 02:07:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (139->64), mult. (177->53), div. (0->0), fcn. (143->6), ass. (0->26)
t17 = cos(qJ(5));
t43 = -m(6) * pkin(4) - m(7) * (t17 * pkin(5) + pkin(4)) - mrSges(5,1);
t42 = m(6) * pkin(8) - m(7) * (-qJ(6) - pkin(8)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t41 = m(7) * pkin(5);
t40 = -m(1) - m(2);
t39 = -mrSges(6,1) - mrSges(7,1);
t38 = mrSges(6,2) + mrSges(7,2);
t37 = -m(5) - m(6) - m(7);
t36 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t15 = sin(qJ(4));
t18 = cos(qJ(4));
t35 = t43 * t15 + t42 * t18 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t14 = sin(qJ(5));
t16 = sin(qJ(1));
t34 = t16 * t14;
t33 = t16 * t17;
t19 = cos(qJ(1));
t32 = t19 * t14;
t31 = t19 * t17;
t12 = pkin(6) + r_base(3);
t30 = pkin(2) + t12;
t29 = t19 * pkin(1) + t16 * qJ(2) + r_base(1);
t27 = t19 * qJ(3) + t29;
t25 = t16 * pkin(1) - t19 * qJ(2) + r_base(2);
t24 = t16 * qJ(3) + t25;
t1 = (-m(1) * r_base(3) - m(4) * t30 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t12 + t37 * (pkin(3) + t30) + (t38 * t14 + t39 * t17 + t43) * t18 - t42 * t15) * g(3) + (-t32 * t41 - m(3) * t25 - m(4) * t24 - mrSges(1,2) + t40 * r_base(2) + t37 * (t19 * pkin(7) + t24) + t39 * (t15 * t33 + t32) - t38 * (-t15 * t34 + t31) - t36 * t19 + t35 * t16) * g(2) + (t34 * t41 - m(3) * t29 - m(4) * t27 - mrSges(1,1) + t40 * r_base(1) + t39 * (t15 * t31 - t34) - t38 * (-t15 * t32 - t33) + t37 * (-t16 * pkin(7) + t27) + t36 * t16 + t35 * t19) * g(1);
U  = t1;
