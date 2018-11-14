% Calculate potential energy for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:22
% EndTime: 2018-11-14 13:45:22
% DurationCPUTime: 0.34s
% Computational Cost: add. (207->60), mult. (218->56), div. (0->0), fcn. (177->10), ass. (0->34)
t51 = -m(1) - m(2);
t50 = -m(4) - m(5);
t20 = sin(pkin(6));
t24 = sin(qJ(1));
t25 = cos(qJ(1));
t38 = pkin(4) - pkin(6);
t32 = cos(t38) / 0.2e1;
t37 = pkin(4) + pkin(6);
t34 = cos(t37);
t30 = t32 + t34 / 0.2e1;
t3 = t24 * t20 - t25 * t30;
t22 = cos(pkin(6));
t31 = sin(t37) / 0.2e1;
t33 = sin(t38);
t29 = t31 - t33 / 0.2e1;
t4 = t24 * t22 + t25 * t29;
t49 = t4 * pkin(2) + t3 * qJ(3);
t48 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3);
t47 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t46 = -m(5) * pkin(3) - t48;
t45 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t21 = sin(pkin(4));
t41 = qJ(2) * t21;
t40 = pkin(5) + r_base(3);
t39 = t24 * pkin(1) + r_base(2);
t23 = cos(pkin(4));
t36 = t23 * qJ(2) + t40;
t35 = t25 * pkin(1) + t24 * t41 + r_base(1);
t28 = -t25 * t41 + t39;
t11 = t32 - t34 / 0.2e1;
t10 = t31 + t33 / 0.2e1;
t6 = t25 * t22 - t24 * t29;
t5 = t25 * t20 + t24 * t30;
t1 = (-m(1) * r_base(3) - m(2) * t40 - m(3) * t36 - mrSges(1,3) - mrSges(2,3) + t50 * (t11 * pkin(2) - t10 * qJ(3) + t36) + t46 * t23 + t45 * t11 - t47 * t10) * g(3) + (-mrSges(1,2) - t24 * mrSges(2,1) - m(3) * t28 - m(4) * (t28 + t49) - m(5) * (t39 + t49) + t51 * r_base(2) + t45 * t4 + t47 * t3 + (-mrSges(2,2) + (-m(5) * (-pkin(3) - qJ(2)) + t48) * t21) * t25) * g(2) + (-m(3) * t35 - t25 * mrSges(2,1) - mrSges(1,1) + t51 * r_base(1) + t50 * (t6 * pkin(2) + t5 * qJ(3) + t35) + t45 * t6 + t47 * t5 + (t46 * t21 + mrSges(2,2)) * t24) * g(1);
U  = t1;
