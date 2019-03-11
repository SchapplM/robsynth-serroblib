% Calculate potential energy for
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:02
% EndTime: 2019-03-09 03:11:02
% DurationCPUTime: 0.48s
% Computational Cost: add. (230->74), mult. (214->68), div. (0->0), fcn. (184->8), ass. (0->33)
t58 = mrSges(4,2) - mrSges(5,3);
t57 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t56 = -m(1) - m(2);
t55 = -m(6) - m(7);
t23 = qJ(1) + pkin(9);
t16 = sin(t23);
t25 = sin(qJ(3));
t43 = qJ(4) * t25;
t28 = cos(qJ(3));
t49 = t16 * t28;
t54 = pkin(3) * t49 + t16 * t43;
t53 = mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t52 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t51 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t50 = t58 * t25 + t57 * t28 - mrSges(3,1);
t17 = cos(t23);
t48 = t17 * t28;
t24 = sin(qJ(5));
t47 = t24 * t25;
t27 = cos(qJ(5));
t46 = t25 * t27;
t42 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t41 = t26 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t40 = t29 * pkin(1) + r_base(1);
t18 = qJ(2) + t42;
t39 = t16 * pkin(2) + t41;
t38 = t25 * pkin(3) + t18;
t37 = t17 * pkin(2) + t16 * pkin(7) + t40;
t33 = -pkin(7) * t17 + t39;
t32 = pkin(3) * t48 + t17 * t43 + t37;
t1 = (-m(1) * r_base(3) - m(2) * t42 - m(5) * t38 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t55 * (t25 * pkin(8) + t38) + (-m(3) - m(4)) * t18 + (-t51 * t27 + t52 * t24 + (m(5) - t55) * qJ(4) - t58) * t28 + t57 * t25) * g(3) + (-mrSges(1,2) - t26 * mrSges(2,1) - t29 * mrSges(2,2) - m(3) * t41 - m(4) * t33 - m(5) * (t33 + t54) + t56 * r_base(2) + t55 * (pkin(8) * t49 + (-pkin(4) - pkin(7)) * t17 + t39 + t54) - t52 * (t16 * t47 - t17 * t27) + t51 * (t16 * t46 + t17 * t24) - t53 * t17 + t50 * t16) * g(2) + (-m(3) * t40 - m(4) * t37 - m(5) * t32 - t29 * mrSges(2,1) + t26 * mrSges(2,2) - mrSges(1,1) + t56 * r_base(1) + t55 * (t16 * pkin(4) + pkin(8) * t48 + t32) - t52 * (t16 * t27 + t17 * t47) - t51 * (t16 * t24 - t17 * t46) + t53 * t16 + t50 * t17) * g(1);
U  = t1;
