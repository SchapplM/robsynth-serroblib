% Calculate potential energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR16_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:39
% EndTime: 2019-12-31 18:38:40
% DurationCPUTime: 0.36s
% Computational Cost: add. (106->62), mult. (142->48), div. (0->0), fcn. (110->6), ass. (0->21)
t36 = -mrSges(4,1) + mrSges(5,2);
t35 = -m(1) - m(2);
t34 = m(5) + m(6);
t11 = sin(qJ(5));
t14 = cos(qJ(5));
t33 = -t11 * mrSges(6,1) - t14 * mrSges(6,2) - mrSges(5,3);
t32 = m(6) * pkin(7) + mrSges(6,3);
t12 = sin(qJ(3));
t15 = cos(qJ(3));
t31 = -t15 * mrSges(4,2) + t36 * t12 + mrSges(2,2) - mrSges(3,3);
t30 = -m(6) * pkin(4) - t14 * mrSges(6,1) + t11 * mrSges(6,2) - mrSges(2,1) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t29 = pkin(3) * t12;
t10 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t27 = t13 * pkin(1) + r_base(2);
t26 = pkin(2) + t10;
t25 = t13 * pkin(6) + t27;
t16 = cos(qJ(1));
t24 = t16 * pkin(1) + t13 * qJ(2) + r_base(1);
t22 = t16 * pkin(6) + t24;
t1 = (-m(1) * r_base(3) - m(4) * t26 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t34 * (t15 * pkin(3) + t12 * qJ(4) + t26) + (-m(2) - m(3)) * t10 + (-t32 + t36) * t15 + (mrSges(4,2) + t33) * t12) * g(3) + (-m(3) * t27 - m(4) * t25 - mrSges(1,2) + t35 * r_base(2) - t34 * (t16 * t15 * qJ(4) + t25) + (m(5) * t29 - (m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3)) * t12 + t33 * t15 + (m(3) + m(4) + t34) * qJ(2) - t31) * t16 + t30 * t13) * g(2) + (-m(3) * t24 - m(4) * t22 - mrSges(1,1) + t35 * r_base(1) - t34 * (t13 * t29 + t22) + t30 * t16 + (-t32 * t12 + (t34 * qJ(4) - t33) * t15 + t31) * t13) * g(1);
U = t1;
