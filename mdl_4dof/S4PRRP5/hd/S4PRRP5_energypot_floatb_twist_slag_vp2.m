% Calculate potential energy for
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:48
% EndTime: 2019-12-31 16:28:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (92->48), mult. (135->45), div. (0->0), fcn. (113->6), ass. (0->20)
t15 = cos(qJ(3));
t34 = -m(4) * pkin(2) - m(5) * (pkin(3) * t15 + pkin(2)) - mrSges(3,1);
t33 = m(4) * pkin(5) - m(5) * (-qJ(4) - pkin(5)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t32 = m(5) * pkin(3);
t31 = -m(1) - m(2);
t30 = -mrSges(4,1) - mrSges(5,1);
t29 = mrSges(2,2) - mrSges(3,3);
t28 = mrSges(4,2) + mrSges(5,2);
t27 = -m(3) - m(4) - m(5);
t14 = sin(qJ(2));
t16 = cos(qJ(2));
t26 = -t33 * t14 + t34 * t16 - mrSges(2,1);
t10 = sin(pkin(6));
t13 = sin(qJ(3));
t25 = t10 * t13;
t11 = cos(pkin(6));
t24 = t11 * t13;
t23 = t13 * t16;
t22 = t15 * t16;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t27) * (qJ(1) + r_base(3)) + t33 * t16 + (t28 * t13 + t30 * t15 + t34) * t14) * g(3) + (t24 * t32 - mrSges(1,2) + t31 * r_base(2) + t27 * (t10 * pkin(1) - t11 * pkin(4) + r_base(2)) + t30 * (t10 * t22 - t24) - t29 * t11 - t28 * (-t10 * t23 - t11 * t15) + t26 * t10) * g(2) + (-t25 * t32 - mrSges(1,1) + t31 * r_base(1) + t30 * (t11 * t22 + t25) - t28 * (t10 * t15 - t11 * t23) + t27 * (t11 * pkin(1) + t10 * pkin(4) + r_base(1)) + t29 * t10 + t26 * t11) * g(1);
U = t1;
