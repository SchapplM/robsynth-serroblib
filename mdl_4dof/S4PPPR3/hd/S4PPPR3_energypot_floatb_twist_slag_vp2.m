% Calculate potential energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:27
% EndTime: 2018-11-14 13:56:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (54->41), mult. (40->24), div. (0->0), fcn. (10->4), ass. (0->11)
t12 = -m(3) - m(4);
t11 = -qJ(3) - pkin(1);
t10 = -m(1) - m(2) - m(5);
t9 = qJ(1) + r_base(2);
t8 = pkin(2) + t9;
t7 = cos(pkin(5));
t6 = sin(pkin(5));
t5 = pkin(5) + qJ(4);
t2 = cos(t5);
t1 = sin(t5);
t3 = (-mrSges(1,3) + mrSges(2,1) + m(3) * pkin(1) - mrSges(3,2) - m(4) * t11 + mrSges(4,3) - m(5) * (-pkin(4) + t11) + mrSges(5,3) + (t10 + t12) * r_base(3)) * g(3) + (-m(1) * r_base(2) - mrSges(1,2) - mrSges(2,3) - mrSges(3,1) - m(4) * t8 - t7 * mrSges(4,1) + t6 * mrSges(4,2) - m(5) * (pkin(3) * t7 + t8) - t2 * mrSges(5,1) + t1 * mrSges(5,2) + (-m(2) - m(3)) * t9) * g(2) + (-mrSges(1,1) + mrSges(2,2) - mrSges(3,3) - t6 * mrSges(4,1) - t7 * mrSges(4,2) - m(5) * (pkin(3) * t6 + qJ(2)) - t1 * mrSges(5,1) - t2 * mrSges(5,2) + t10 * r_base(1) + t12 * (qJ(2) + r_base(1))) * g(1);
U  = t3;
