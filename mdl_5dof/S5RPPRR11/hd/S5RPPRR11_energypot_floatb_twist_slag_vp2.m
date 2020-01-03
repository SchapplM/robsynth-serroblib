% Calculate potential energy for
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:29
% DurationCPUTime: 0.31s
% Computational Cost: add. (98->55), mult. (121->38), div. (0->0), fcn. (89->6), ass. (0->21)
t11 = cos(qJ(5));
t8 = sin(qJ(5));
t31 = -m(6) * pkin(4) - t11 * mrSges(6,1) + t8 * mrSges(6,2) - mrSges(5,1);
t30 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t29 = -m(1) - m(2);
t28 = m(5) + m(6);
t12 = cos(qJ(4));
t9 = sin(qJ(4));
t26 = -t30 * t12 + t31 * t9 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t25 = t8 * mrSges(6,1) + t11 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t7 = pkin(5) + r_base(3);
t10 = sin(qJ(1));
t24 = t10 * pkin(1) + r_base(2);
t23 = pkin(2) + t7;
t13 = cos(qJ(1));
t22 = t13 * pkin(1) + t10 * qJ(2) + r_base(1);
t17 = -qJ(2) * t13 + t24;
t1 = t10 * qJ(3);
t16 = t1 + t17;
t5 = t13 * pkin(6);
t2 = (-m(1) * r_base(3) - m(4) * t23 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - t28 * (pkin(3) + t23) + t30 * t9 + (-m(2) - m(3)) * t7 + t31 * t12) * g(3) + (-mrSges(1,2) - m(3) * t17 - m(4) * t16 - m(5) * (t16 + t5) - m(6) * (t1 + t5 + t24) + t29 * r_base(2) + (m(6) * qJ(2) - t25) * t13 + t26 * t10) * g(2) + (-m(3) * t22 - mrSges(1,1) + t29 * r_base(1) + (-m(4) - t28) * (t13 * qJ(3) + t22) + t26 * t13 + (t28 * pkin(6) + t25) * t10) * g(1);
U = t2;
