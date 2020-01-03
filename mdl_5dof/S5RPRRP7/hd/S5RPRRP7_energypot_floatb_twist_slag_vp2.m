% Calculate potential energy for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:29
% EndTime: 2019-12-31 18:44:30
% DurationCPUTime: 0.41s
% Computational Cost: add. (173->53), mult. (166->41), div. (0->0), fcn. (142->8), ass. (0->23)
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t41 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t42 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t53 = t41 * t20 + t42 * t23 - mrSges(4,1);
t50 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t45 = -m(5) - m(6);
t49 = -m(4) + t45;
t48 = t42 * t20 - t41 * t23 + mrSges(3,2) - mrSges(4,3);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t47 = -mrSges(3,1) + (t45 * pkin(7) + t50) * t21 + (t45 * pkin(3) + t53) * t24;
t46 = -m(1) - m(2);
t34 = pkin(5) + r_base(3);
t22 = sin(qJ(1));
t33 = t22 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t32 = t25 * pkin(1) + r_base(1);
t15 = qJ(2) + t34;
t19 = qJ(1) + pkin(8);
t14 = cos(t19);
t13 = sin(t19);
t1 = (-m(1) * r_base(3) - m(2) * t34 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t45 * (t21 * pkin(3) - t24 * pkin(7) + t15) + (-m(3) - m(4)) * t15 - t50 * t24 + t53 * t21) * g(3) + (-m(3) * t33 - t22 * mrSges(2,1) - t25 * mrSges(2,2) + t46 * r_base(2) - mrSges(1,2) + t49 * (t13 * pkin(2) - t14 * pkin(6) + t33) - t48 * t14 + t47 * t13) * g(2) + (-m(3) * t32 - t25 * mrSges(2,1) + t22 * mrSges(2,2) + t46 * r_base(1) - mrSges(1,1) + t49 * (t14 * pkin(2) + t13 * pkin(6) + t32) + t48 * t13 + t47 * t14) * g(1);
U = t1;
