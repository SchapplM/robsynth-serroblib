% Calculate potential energy for
% S2PP1
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
%   pkin=[a2]';
% m [3x1]
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
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S2PP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2PP1_energypot_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-03-03 18:41:17
% EndTime: 2021-03-03 18:41:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (21->17), mult. (18->9), div. (0->0), fcn. (0->0), ass. (0->3)
t3 = -m(2) - m(3);
t2 = -m(1) + t3;
t1 = (-m(3) * pkin(1) + t2 * r_base(3) - mrSges(2,1) + mrSges(3,2) - mrSges(1,3)) * g(3) + (-m(3) * qJ(2) + t2 * r_base(2) - mrSges(1,2) + mrSges(2,2) - mrSges(3,3)) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(3,1) - mrSges(2,3) + t3 * (qJ(1) + r_base(1))) * g(1);
U = t1;
