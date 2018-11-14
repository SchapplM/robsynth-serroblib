% Calculate potential energy for
% S4RRPP1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RRPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:28
% EndTime: 2018-11-14 13:51:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (92->41), mult. (60->31), div. (0->0), fcn. (30->6), ass. (0->17)
t21 = -m(4) - m(5);
t11 = qJ(1) + qJ(2);
t18 = pkin(4) + r_base(3);
t17 = pkin(5) + t18;
t16 = -m(1) - m(2) - m(3) + t21;
t15 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t14 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t13 = cos(qJ(1));
t12 = sin(qJ(1));
t10 = t13 * pkin(1);
t9 = t12 * pkin(1);
t8 = cos(t11);
t7 = sin(t11);
t6 = pkin(6) + t11;
t2 = cos(t6);
t1 = sin(t6);
t3 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t18 - mrSges(2,3) - m(3) * t17 - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) + t21 * (qJ(3) + t17)) * g(3) + (-m(3) * t9 - t12 * mrSges(2,1) - t7 * mrSges(3,1) - mrSges(2,2) * t13 - t8 * mrSges(3,2) + t15 * t1 + t14 * t2 + t16 * r_base(2) - mrSges(1,2) + t21 * (pkin(2) * t7 + t9)) * g(2) + (-m(3) * t10 - mrSges(2,1) * t13 - t8 * mrSges(3,1) + t12 * mrSges(2,2) + t7 * mrSges(3,2) - t14 * t1 + t15 * t2 + t16 * r_base(1) - mrSges(1,1) + t21 * (pkin(2) * t8 + t10)) * g(1);
U  = t3;
