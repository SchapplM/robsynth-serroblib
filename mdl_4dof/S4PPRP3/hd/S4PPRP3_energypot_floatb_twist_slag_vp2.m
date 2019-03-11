% Calculate potential energy for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:29
% EndTime: 2019-03-08 18:14:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (50->37), mult. (40->18), div. (0->0), fcn. (10->2), ass. (0->9)
t11 = -m(3) - m(4);
t10 = -pkin(4) - pkin(1);
t9 = mrSges(4,2) + mrSges(5,2);
t8 = -m(1) - m(2) - m(5);
t7 = qJ(1) + r_base(2);
t5 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t4 = cos(qJ(3));
t3 = sin(qJ(3));
t1 = (-mrSges(1,3) + mrSges(2,1) + m(3) * pkin(1) - mrSges(3,2) - m(4) * t10 + mrSges(4,3) - m(5) * (-qJ(4) + t10) + mrSges(5,3) + (t8 + t11) * r_base(3)) * g(3) + (-m(1) * r_base(2) - mrSges(1,2) - mrSges(2,3) - mrSges(3,1) + t5 * t4 + t9 * t3 + (-m(2) - m(3)) * t7 + (-m(4) - m(5)) * (pkin(2) + t7)) * g(2) + (-m(5) * qJ(2) - mrSges(1,1) + mrSges(2,2) - mrSges(3,3) + t8 * r_base(1) - t9 * t4 + t5 * t3 + t11 * (qJ(2) + r_base(1))) * g(1);
U  = t1;
