% Calculate potential energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3RRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:47
% EndTime: 2018-11-14 10:15:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (55->35), mult. (42->29), div. (0->0), fcn. (18->6), ass. (0->14)
t8 = qJ(1) + qJ(2);
t13 = pkin(3) + r_base(3);
t12 = -m(1) - m(2) - m(3) - m(4);
t11 = pkin(4) + t13;
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t7 = t10 * pkin(1);
t6 = t9 * pkin(1);
t5 = qJ(3) + t8;
t4 = cos(t8);
t3 = sin(t8);
t2 = cos(t5);
t1 = sin(t5);
t14 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t13 - mrSges(2,3) - m(3) * t11 - mrSges(3,3) - m(4) * (pkin(5) + t11) - mrSges(4,3)) * g(3) + (-mrSges(1,2) - t9 * mrSges(2,1) - t10 * mrSges(2,2) - m(3) * t6 - t3 * mrSges(3,1) - t4 * mrSges(3,2) - m(4) * (pkin(2) * t3 + t6) - t1 * mrSges(4,1) - t2 * mrSges(4,2) + t12 * r_base(2)) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t10 + t9 * mrSges(2,2) - m(3) * t7 - t4 * mrSges(3,1) + t3 * mrSges(3,2) - m(4) * (pkin(2) * t4 + t7) - t2 * mrSges(4,1) + t1 * mrSges(4,2) + t12 * r_base(1)) * g(1);
U  = t14;
