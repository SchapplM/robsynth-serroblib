% Calculate potential energy for
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:28
% EndTime: 2019-03-09 16:43:29
% DurationCPUTime: 0.48s
% Computational Cost: add. (222->77), mult. (228->68), div. (0->0), fcn. (198->8), ass. (0->36)
t59 = -mrSges(4,1) + mrSges(5,2);
t58 = mrSges(4,2) - mrSges(5,3);
t57 = -m(2) - m(3);
t56 = -m(6) - m(7);
t55 = -mrSges(6,3) - mrSges(7,2);
t54 = -m(1) + t57;
t23 = qJ(2) + qJ(3);
t18 = sin(t23);
t19 = cos(t23);
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t53 = -m(3) * pkin(1) - t28 * mrSges(3,1) + t25 * mrSges(3,2) + t58 * t18 + t59 * t19 - mrSges(2,1);
t52 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t51 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t50 = m(3) * pkin(7) + mrSges(5,1) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t26 = sin(qJ(1));
t49 = t19 * t26;
t29 = cos(qJ(1));
t48 = t19 * t29;
t24 = sin(qJ(5));
t47 = t24 * t26;
t46 = t24 * t29;
t27 = cos(qJ(5));
t45 = t26 * t27;
t44 = t27 * t29;
t43 = qJ(4) * t18;
t22 = pkin(6) + r_base(3);
t42 = t25 * pkin(2) + t22;
t16 = pkin(2) * t28 + pkin(1);
t30 = -pkin(8) - pkin(7);
t41 = t26 * t16 + t29 * t30 + r_base(2);
t40 = t18 * pkin(3) + t42;
t37 = t29 * t16 - t26 * t30 + r_base(1);
t36 = pkin(3) * t49 + t26 * t43 + t41;
t34 = pkin(3) * t48 + t29 * t43 + t37;
t1 = (-m(1) * r_base(3) - m(4) * t42 - m(5) * t40 - t25 * mrSges(3,1) - t28 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t56 * (t18 * pkin(9) + t40) + t57 * t22 + (-t51 * t27 + t52 * t24 + (m(5) - t56) * qJ(4) - t58) * t19 + (t55 + t59) * t18) * g(3) + (-m(4) * t41 - m(5) * t36 - mrSges(1,2) + t55 * t49 + t56 * (-pkin(4) * t29 + pkin(9) * t49 + t36) - t52 * (t18 * t47 - t44) + t51 * (t18 * t45 + t46) + t54 * r_base(2) + t50 * t29 + t53 * t26) * g(2) + (-m(4) * t37 - m(5) * t34 - mrSges(1,1) + t55 * t48 + t56 * (t26 * pkin(4) + pkin(9) * t48 + t34) - t52 * (t18 * t46 + t45) - t51 * (-t18 * t44 + t47) + t54 * r_base(1) + t53 * t29 - t50 * t26) * g(1);
U  = t1;
