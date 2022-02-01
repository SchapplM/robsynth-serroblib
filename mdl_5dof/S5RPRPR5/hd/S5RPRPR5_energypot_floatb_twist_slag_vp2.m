% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:47
% EndTime: 2022-01-23 09:24:48
% DurationCPUTime: 0.58s
% Computational Cost: add. (170->73), mult. (185->63), div. (0->0), fcn. (161->10), ass. (0->31)
t24 = cos(qJ(3));
t43 = t24 * pkin(3);
t17 = qJ(3) + pkin(9);
t10 = qJ(5) + t17;
t5 = sin(t10);
t6 = cos(t10);
t8 = sin(t17);
t9 = cos(t17);
t59 = -t9 * mrSges(5,1) + t8 * mrSges(5,2) - m(5) * t43 - mrSges(3,1) - m(6) * (pkin(4) * t9 + pkin(2) + t43) - t6 * mrSges(6,1) + t5 * mrSges(6,2);
t21 = qJ(4) + pkin(6);
t58 = -mrSges(4,3) - mrSges(5,3) - m(4) * pkin(6) - m(5) * t21 + mrSges(3,2) + m(6) * (-pkin(7) - t21) - mrSges(6,3);
t57 = -m(4) - m(5);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t54 = t57 * pkin(1) - mrSges(2,1) + (t57 * pkin(2) + t59) * t20 + t58 * t19;
t53 = -m(3) - m(6);
t51 = -m(1) - m(2) - m(5);
t22 = sin(qJ(3));
t44 = t22 * pkin(3);
t45 = m(5) * (qJ(2) + t44) + t8 * mrSges(5,1) + t9 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + t5 * mrSges(6,1) + t6 * mrSges(6,2);
t23 = sin(qJ(1));
t38 = t23 * t22;
t37 = t23 * t24;
t25 = cos(qJ(1));
t35 = t25 * t22;
t34 = t25 * t24;
t18 = pkin(5) + r_base(3);
t30 = -t25 * qJ(2) + r_base(2);
t14 = t23 * pkin(1);
t3 = pkin(4) * t8 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t57 * (t19 * pkin(2) + t18) + (-m(2) + t53) * t18 - t58 * t20 + (-t24 * mrSges(4,1) + t22 * mrSges(4,2) + t59) * t19) * g(3) + (-mrSges(1,2) - m(3) * (t14 + t30) - m(4) * t30 - (t20 * t37 - t35) * mrSges(4,1) - (-t20 * t38 - t34) * mrSges(4,2) - m(6) * t14 + (-m(6) + t51) * r_base(2)) * g(2) + (-mrSges(1,1) - (t20 * t34 + t38) * mrSges(4,1) - (-t20 * t35 + t37) * mrSges(4,2) + t51 * r_base(1) + (-m(4) + t53) * (t23 * qJ(2) + r_base(1))) * g(1) + ((-m(6) * (-qJ(2) - t3) + t45) * g(2) + (t53 * pkin(1) + t54) * g(1)) * t25 + (t54 * g(2) + (-m(6) * t3 - t45) * g(1)) * t23;
U = t1;
