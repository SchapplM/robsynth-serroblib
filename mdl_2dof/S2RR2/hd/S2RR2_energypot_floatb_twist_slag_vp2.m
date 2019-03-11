% Calculate potential energy for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S2RR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR2_energypot_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:52
% EndTime: 2019-03-08 18:00:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (26->17), mult. (36->15), div. (0->0), fcn. (18->4), ass. (0->8)
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t8 = t3 * mrSges(3,1) - t1 * mrSges(3,2) + mrSges(2,1);
t7 = -m(1) - m(2) - m(3);
t6 = -m(3) * pkin(1) + mrSges(2,2) - mrSges(3,3);
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t5 = (t8 * t2 + t6 * t4 + t7 * r_base(3) - mrSges(1,3)) * g(3) + (-mrSges(3,1) * t1 - mrSges(3,2) * t3 + t7 * r_base(2) - mrSges(1,2) - mrSges(2,3)) * g(2) + (t6 * t2 - t8 * t4 + t7 * r_base(1) - mrSges(1,1)) * g(1);
U  = t5;
