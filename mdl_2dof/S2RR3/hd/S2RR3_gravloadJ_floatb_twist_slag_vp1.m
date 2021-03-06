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
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2RR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [3x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:24
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (21->12), mult. (23->19), div. (0->0), fcn. (14->4), ass. (0->8)
t4 = qJ(1) + qJ(2);
t2 = sin(t4);
t3 = cos(t4);
t8 = t3 * rSges(3,1) - t2 * rSges(3,2);
t7 = -t2 * rSges(3,1) - t3 * rSges(3,2);
t6 = cos(qJ(1));
t5 = sin(qJ(1));
t1 = [-m(2) * (g(1) * (-t5 * rSges(2,1) - t6 * rSges(2,2)) + g(2) * (t6 * rSges(2,1) - t5 * rSges(2,2))) - m(3) * (g(1) * (-t5 * pkin(1) + t7) + g(2) * (t6 * pkin(1) + t8)), -m(3) * (g(1) * t7 + g(2) * t8)];
taug = t1(:);
