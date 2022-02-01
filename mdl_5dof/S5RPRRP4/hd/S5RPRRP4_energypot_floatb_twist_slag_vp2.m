% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:42
% EndTime: 2022-01-23 09:31:43
% DurationCPUTime: 0.55s
% Computational Cost: add. (160->70), mult. (185->63), div. (0->0), fcn. (161->8), ass. (0->32)
t26 = pkin(7) + pkin(6);
t60 = -m(4) * pkin(6) - m(5) * t26 + m(6) * (-qJ(5) - t26) + mrSges(3,2) - mrSges(6,3) - mrSges(4,3) - mrSges(5,3);
t59 = -m(4) - m(5);
t19 = qJ(3) + qJ(4);
t11 = cos(t19);
t24 = cos(qJ(3));
t45 = t24 * pkin(3);
t58 = -m(5) * t45 - m(6) * (pkin(4) * t11 + pkin(2) + t45) - mrSges(3,1);
t54 = -m(3) - m(6);
t56 = -m(4) + t54;
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t55 = t54 * pkin(1) - mrSges(2,1) + t59 * (pkin(2) * t21 + pkin(1)) + t58 * t21 + t60 * t20;
t53 = -mrSges(5,1) - mrSges(6,1);
t52 = mrSges(5,2) + mrSges(6,2);
t50 = -m(1) - m(2) - m(5);
t10 = sin(t19);
t22 = sin(qJ(3));
t46 = t22 * pkin(3);
t48 = m(5) * (qJ(2) + t46) + m(6) * (pkin(4) * t10 + t46) - mrSges(2,2) + mrSges(3,3);
t23 = sin(qJ(1));
t43 = t23 * t10;
t42 = t23 * t11;
t40 = t23 * t22;
t39 = t23 * t24;
t25 = cos(qJ(1));
t38 = t25 * t10;
t37 = t25 * t11;
t35 = t25 * t22;
t34 = t25 * t24;
t18 = pkin(5) + r_base(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t59 * (t20 * pkin(2) + t18) + (-m(2) + t54) * t18 - t60 * t21 + (-t24 * mrSges(4,1) + t22 * mrSges(4,2) + t52 * t10 + t53 * t11 + t58) * t20) * g(3) + (-mrSges(1,2) - (t21 * t39 - t35) * mrSges(4,1) - (-t21 * t40 - t34) * mrSges(4,2) + t53 * (t21 * t42 - t38) + t50 * r_base(2) + t48 * t25 - t52 * (-t21 * t43 - t37) + t56 * (-t25 * qJ(2) + r_base(2)) + t55 * t23) * g(2) + (-mrSges(1,1) - (t21 * t34 + t40) * mrSges(4,1) - (-t21 * t35 + t39) * mrSges(4,2) + t53 * (t21 * t37 + t43) - t52 * (-t21 * t38 + t42) + t50 * r_base(1) - t48 * t23 + t56 * (t23 * qJ(2) + r_base(1)) + t55 * t25) * g(1);
U = t1;
