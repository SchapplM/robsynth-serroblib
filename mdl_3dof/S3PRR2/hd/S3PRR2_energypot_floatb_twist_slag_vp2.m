% Calculate potential energy for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3PRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:44
% EndTime: 2018-11-14 10:12:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (41->31), mult. (34->22), div. (0->0), fcn. (10->4), ass. (0->9)
t9 = -m(2) - m(3);
t8 = qJ(1) + r_base(1);
t7 = -m(1) - m(4) + t9;
t6 = cos(qJ(2));
t5 = sin(qJ(2));
t4 = qJ(2) + qJ(3);
t2 = cos(t4);
t1 = sin(t4);
t3 = (-mrSges(1,3) - mrSges(2,2) + m(3) * pkin(3) + mrSges(3,3) - m(4) * (-pkin(4) - pkin(3)) + mrSges(4,3) + t7 * r_base(3)) * g(3) + (-mrSges(1,2) - mrSges(2,1) - m(3) * pkin(1) - t6 * mrSges(3,1) + t5 * mrSges(3,2) - m(4) * (pkin(2) * t6 + pkin(1)) - t2 * mrSges(4,1) + t1 * mrSges(4,2) + t7 * r_base(2)) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) - t5 * mrSges(3,1) - t6 * mrSges(3,2) - m(4) * (pkin(2) * t5 + t8) - t1 * mrSges(4,1) - t2 * mrSges(4,2) + t9 * t8) * g(1);
U  = t3;
