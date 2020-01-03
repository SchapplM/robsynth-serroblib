% Calculate potential energy for
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:19
% EndTime: 2019-12-31 16:57:20
% DurationCPUTime: 0.23s
% Computational Cost: add. (96->43), mult. (98->32), div. (0->0), fcn. (68->6), ass. (0->18)
t27 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t26 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t25 = -m(2) - m(3);
t24 = -m(4) - m(5);
t23 = -m(1) + t25;
t11 = sin(qJ(2));
t13 = cos(qJ(2));
t8 = qJ(2) + pkin(6);
t5 = sin(t8);
t6 = cos(t8);
t22 = -m(3) * pkin(1) - t13 * mrSges(3,1) + t11 * mrSges(3,2) - t26 * t5 + t27 * t6 - mrSges(2,1);
t21 = m(3) * pkin(5) - mrSges(2,2) + mrSges(5,2) + mrSges(3,3) + mrSges(4,3);
t9 = pkin(4) + r_base(3);
t14 = cos(qJ(1));
t12 = sin(qJ(1));
t10 = -qJ(3) - pkin(5);
t4 = pkin(2) * t13 + pkin(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) * t11 - mrSges(3,2) * t13 - mrSges(1,3) - mrSges(2,3) + t25 * t9 + t24 * (t11 * pkin(2) + t9) + t26 * t6 + t27 * t5) * g(3) + (-mrSges(1,2) + t24 * (t14 * t10 + t12 * t4 + r_base(2)) + t23 * r_base(2) + t21 * t14 + t22 * t12) * g(2) + (-mrSges(1,1) + t24 * (-t12 * t10 + t14 * t4 + r_base(1)) + t23 * r_base(1) + t22 * t14 - t21 * t12) * g(1);
U = t1;
