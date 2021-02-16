% Calculate potential energy for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [2x1]
%   mass of all robot links (including the base)
% mrSges [2x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S1P1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(2,1),zeros(2,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1P1_energypot_floatb_twist_slag_vp2: mrSges has to be [2x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:50
% EndTime: 2021-01-14 12:21:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (12->11), mult. (12->7), div. (0->0), fcn. (0->0), ass. (0->2)
t1 = -m(1) - m(2);
t2 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * (qJ(1) + r_base(3)) - mrSges(2,3)) * g(3) + (t1 * r_base(2) - mrSges(1,2) - mrSges(2,2)) * g(2) + (t1 * r_base(1) - mrSges(1,1) - mrSges(2,1)) * g(1);
U = t2;
