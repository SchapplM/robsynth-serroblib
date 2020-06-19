% Calculate Gravitation load on the joints for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2RR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [3x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:24
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (20->10), mult. (22->13), div. (0->0), fcn. (14->4), ass. (0->9)
t9 = m(3) * pkin(1) + mrSges(2,1);
t4 = qJ(1) + qJ(2);
t2 = sin(t4);
t3 = cos(t4);
t8 = mrSges(3,1) * t3 - t2 * mrSges(3,2);
t7 = -t2 * mrSges(3,1) - t3 * mrSges(3,2);
t6 = cos(qJ(1));
t5 = sin(qJ(1));
t1 = [(t5 * mrSges(2,2) - t6 * t9 - t8) * g(2) + (t6 * mrSges(2,2) + t5 * t9 - t7) * g(1), -g(1) * t7 - g(2) * t8];
taug = t1(:);
