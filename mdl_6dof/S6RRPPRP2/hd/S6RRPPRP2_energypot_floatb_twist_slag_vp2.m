% Calculate potential energy for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:33
% EndTime: 2019-03-09 08:29:34
% DurationCPUTime: 0.50s
% Computational Cost: add. (217->79), mult. (218->68), div. (0->0), fcn. (184->8), ass. (0->36)
t22 = sin(qJ(5));
t58 = -m(7) * pkin(5) * t22 + mrSges(4,2) - mrSges(5,3);
t57 = m(7) * (-qJ(6) - pkin(8)) - mrSges(4,1) + mrSges(5,2) - mrSges(7,3);
t56 = -m(2) - m(3);
t55 = -m(5) - m(7);
t27 = cos(qJ(1));
t18 = qJ(2) + pkin(9);
t15 = sin(t18);
t39 = qJ(4) * t15;
t16 = cos(t18);
t44 = t16 * t27;
t54 = pkin(3) * t44 + t27 * t39;
t53 = mrSges(6,1) + mrSges(7,1);
t52 = mrSges(6,2) + mrSges(7,2);
t51 = -m(1) + t56;
t50 = m(6) - t55;
t49 = -m(6) * pkin(8) - mrSges(6,3);
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t48 = -m(3) * pkin(1) - t26 * mrSges(3,1) + t23 * mrSges(3,2) + t58 * t15 + t57 * t16 - mrSges(2,1);
t25 = cos(qJ(5));
t47 = -m(3) * pkin(7) - m(7) * (pkin(5) * t25 + pkin(4)) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t24 = sin(qJ(1));
t45 = t16 * t24;
t43 = t22 * t27;
t42 = t24 * t22;
t41 = t24 * t25;
t40 = t25 * t27;
t19 = pkin(6) + r_base(3);
t14 = pkin(2) * t26 + pkin(1);
t38 = t27 * t14 + r_base(1);
t37 = t23 * pkin(2) + t19;
t21 = -qJ(3) - pkin(7);
t36 = t24 * t14 + t27 * t21 + r_base(2);
t31 = -t24 * t21 + t38;
t1 = (-m(1) * r_base(3) - m(4) * t37 - t23 * mrSges(3,1) - t26 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t56 * t19 - t50 * (t15 * pkin(3) + t37) + (qJ(4) * t50 + t22 * t53 + t25 * t52 - t58) * t16 + (t49 + t57) * t15) * g(3) + (-m(4) * t36 - mrSges(1,2) + t49 * t45 - t53 * (t15 * t42 - t40) - t50 * (pkin(3) * t45 + t24 * t39 + t36) - t52 * (t15 * t41 + t43) + t51 * r_base(2) + (m(6) * pkin(4) - t47) * t27 + t48 * t24) * g(2) + (-mrSges(1,1) - m(4) * t31 - m(6) * (pkin(8) * t44 + t38 + t54) - mrSges(6,3) * t44 + t55 * (t31 + t54) - t53 * (t15 * t43 + t41) - t52 * (t15 * t40 - t42) + t51 * r_base(1) + (-m(6) * (pkin(4) - t21) + t47) * t24 + t48 * t27) * g(1);
U  = t1;
