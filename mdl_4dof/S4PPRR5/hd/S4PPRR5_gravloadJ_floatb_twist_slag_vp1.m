% Calculate Gravitation load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (38->22), mult. (86->37), div. (0->0), fcn. (71->6), ass. (0->13)
t2 = sin(pkin(6));
t3 = cos(pkin(6));
t22 = -g(1) * t2 + g(2) * t3;
t4 = sin(qJ(4));
t6 = cos(qJ(4));
t8 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t6 - rSges(5,2) * t4 + pkin(3));
t19 = -m(4) * rSges(4,2) + m(5) * (rSges(5,3) + pkin(5));
t5 = sin(qJ(3));
t16 = t4 * t5;
t15 = t5 * t6;
t11 = -m(3) - m(4) - m(5);
t7 = cos(qJ(3));
t1 = [(-m(2) + t11) * g(3), -t11 * t22, (-t19 * t7 + t8 * t5) * g(3) + t22 * (t19 * t5 + t8 * t7), -m(5) * (g(1) * ((-t2 * t16 + t3 * t6) * rSges(5,1) + (-t2 * t15 - t3 * t4) * rSges(5,2)) + g(2) * ((t3 * t16 + t2 * t6) * rSges(5,1) + (t3 * t15 - t2 * t4) * rSges(5,2)) + g(3) * (-rSges(5,1) * t4 - rSges(5,2) * t6) * t7)];
taug = t1(:);
