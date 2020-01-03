% Calculate potential energy for
% S5RRPRP11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:17
% EndTime: 2019-12-31 20:12:17
% DurationCPUTime: 0.45s
% Computational Cost: add. (126->66), mult. (194->59), div. (0->0), fcn. (170->6), ass. (0->30)
t51 = -mrSges(3,1) + mrSges(4,2);
t50 = mrSges(3,2) - mrSges(4,3);
t49 = -m(1) - m(2);
t48 = -m(5) - m(6);
t21 = sin(qJ(1));
t20 = sin(qJ(2));
t35 = qJ(3) * t20;
t23 = cos(qJ(2));
t39 = t21 * t23;
t47 = pkin(2) * t39 + t21 * t35;
t46 = -mrSges(5,3) - mrSges(6,2);
t45 = t50 * t20 + t51 * t23 - mrSges(2,1);
t44 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t43 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t42 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t19 = sin(qJ(4));
t41 = t21 * t19;
t22 = cos(qJ(4));
t40 = t21 * t22;
t24 = cos(qJ(1));
t38 = t24 * t19;
t37 = t24 * t22;
t36 = t24 * t23;
t18 = pkin(5) + r_base(3);
t34 = t21 * pkin(1) + r_base(2);
t33 = t20 * pkin(2) + t18;
t32 = t24 * pkin(1) + t21 * pkin(6) + r_base(1);
t28 = -t24 * pkin(6) + t34;
t27 = pkin(2) * t36 + t24 * t35 + t32;
t1 = (-m(1) * r_base(3) - m(4) * t33 - mrSges(1,3) - mrSges(2,3) + t48 * (t20 * pkin(7) + t33) + (-m(2) - m(3)) * t18 + (-t42 * t22 + t43 * t19 + (m(4) - t48) * qJ(3) - t50) * t23 + (t46 + t51) * t20) * g(3) + (-mrSges(1,2) - m(3) * t28 - m(4) * (t28 + t47) + t49 * r_base(2) + t48 * (pkin(7) * t39 + (-pkin(3) - pkin(6)) * t24 + t34 + t47) - t43 * (t20 * t41 - t37) + t46 * t39 + t42 * (t20 * t40 + t38) - t44 * t24 + t45 * t21) * g(2) + (-m(3) * t32 - m(4) * t27 - mrSges(1,1) + t49 * r_base(1) + t46 * t36 + t48 * (t21 * pkin(3) + pkin(7) * t36 + t27) - t43 * (t20 * t38 + t40) - t42 * (-t20 * t37 + t41) + t45 * t24 + t44 * t21) * g(1);
U = t1;
