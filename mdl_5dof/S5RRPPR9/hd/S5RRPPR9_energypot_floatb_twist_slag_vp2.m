% Calculate potential energy for
% S5RRPPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:36
% EndTime: 2019-12-31 19:40:36
% DurationCPUTime: 0.37s
% Computational Cost: add. (117->64), mult. (171->51), div. (0->0), fcn. (139->6), ass. (0->27)
t45 = -m(6) * pkin(7) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t44 = t17 * mrSges(6,1) - t14 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t43 = -m(1) - m(2);
t42 = m(5) + m(6);
t16 = sin(qJ(1));
t15 = sin(qJ(2));
t34 = qJ(3) * t15;
t18 = cos(qJ(2));
t36 = t16 * t18;
t41 = pkin(2) * t36 + t16 * t34;
t19 = cos(qJ(1));
t40 = pkin(3) * t36 + t19 * qJ(4);
t38 = -mrSges(2,1) + t45 * t18 + (-m(6) * pkin(4) - t44) * t15;
t37 = t14 * mrSges(6,1) + t17 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t35 = t19 * t18;
t13 = pkin(5) + r_base(3);
t33 = t16 * pkin(1) + r_base(2);
t32 = t15 * pkin(2) + t13;
t31 = t19 * pkin(1) + t16 * pkin(6) + r_base(1);
t25 = -t19 * pkin(6) + t33;
t24 = pkin(2) * t35 + t19 * t34 + t31;
t22 = -t18 * qJ(3) + t32;
t21 = t25 + t41;
t8 = t15 * pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t22 - m(5) * (t22 + t8) - m(6) * (t32 + t8) + (-m(2) - m(3)) * t13 + (-m(6) * (-pkin(4) - qJ(3)) + t44) * t18 + t45 * t15) * g(3) + (-mrSges(1,2) - m(3) * t25 - m(4) * t21 - m(5) * (t21 + t40) - m(6) * (t33 + t40 + t41) + t43 * r_base(2) + (m(6) * pkin(6) - t37) * t19 + t38 * t16) * g(2) + (-m(3) * t31 - m(4) * t24 - mrSges(1,1) + t43 * r_base(1) - t42 * (pkin(3) * t35 + t24) + t38 * t19 + (t42 * qJ(4) + t37) * t16) * g(1);
U = t1;
