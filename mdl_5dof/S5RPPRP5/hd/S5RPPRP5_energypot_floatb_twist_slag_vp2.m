% Calculate potential energy for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.42s
% Computational Cost: add. (130->57), mult. (206->54), div. (0->0), fcn. (188->6), ass. (0->27)
t52 = -mrSges(3,1) - mrSges(4,1);
t51 = mrSges(3,2) - mrSges(4,3);
t24 = sin(qJ(4));
t44 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t50 = t44 * t24;
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t26 = cos(qJ(4));
t41 = t22 * t26;
t45 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t5 = t22 * t24 + t23 * t26;
t49 = t51 * t22 + t52 * t23 - t44 * t41 + t5 * t45 - mrSges(2,1);
t48 = -m(1) - m(2);
t47 = -m(5) - m(6);
t43 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3) + mrSges(6,2);
t25 = sin(qJ(1));
t40 = t25 * t23;
t27 = cos(qJ(1));
t39 = t27 * t23;
t38 = qJ(3) * t22;
t21 = pkin(5) + r_base(3);
t37 = t27 * pkin(1) + t25 * qJ(2) + r_base(1);
t34 = t25 * pkin(1) - t27 * qJ(2) + r_base(2);
t33 = pkin(2) * t39 + t27 * t38 + t37;
t32 = t22 * pkin(2) - t23 * qJ(3) + t21;
t30 = pkin(2) * t40 + t25 * t38 + t34;
t1 = (-m(1) * r_base(3) - m(4) * t32 - mrSges(1,3) - mrSges(2,3) + t47 * (t22 * pkin(3) + t32) + t45 * (-t23 * t24 + t41) + t44 * t5 - t51 * t23 + t52 * t22 + (-m(2) - m(3)) * t21) * g(3) + (-m(3) * t34 - m(4) * t30 - mrSges(1,2) + t48 * r_base(2) + t47 * (pkin(3) * t40 + t27 * pkin(6) + t30) + t40 * t50 - t43 * t27 + t49 * t25) * g(2) + (-m(3) * t37 - m(4) * t33 - mrSges(1,1) + t48 * r_base(1) + t47 * (pkin(3) * t39 - t25 * pkin(6) + t33) + t39 * t50 + t43 * t25 + t49 * t27) * g(1);
U = t1;
