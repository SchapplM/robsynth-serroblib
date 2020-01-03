% Calculate Gravitation load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (108->56), mult. (242->82), div. (0->0), fcn. (243->8), ass. (0->26)
t43 = rSges(6,3) + pkin(6);
t21 = sin(qJ(2));
t23 = cos(qJ(2));
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t6 = t21 * t33 + t23 * t34;
t18 = sin(pkin(7));
t19 = cos(pkin(7));
t42 = g(1) * t19 + g(2) * t18;
t41 = t42 * t21;
t40 = -m(5) - m(6);
t36 = t23 * pkin(2) + t21 * qJ(3);
t35 = qJ(3) * t23;
t32 = -m(4) + t40;
t31 = t23 * pkin(3) + t36;
t25 = -t21 * t34 + t23 * t33;
t24 = (-pkin(2) - pkin(3)) * t41;
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t14 = t19 * t35;
t13 = t18 * t35;
t5 = t25 * t19;
t4 = t6 * t19;
t3 = t25 * t18;
t2 = t6 * t18;
t1 = [(-m(2) - m(3) + t32) * g(3), -m(3) * (g(3) * (rSges(3,1) * t23 - t21 * rSges(3,2)) + t42 * (-rSges(3,1) * t21 - rSges(3,2) * t23)) - m(4) * (g(1) * t14 + g(2) * t13 + g(3) * (rSges(4,1) * t23 + t21 * rSges(4,3) + t36) + t42 * (rSges(4,3) * t23 + (-rSges(4,1) - pkin(2)) * t21)) - m(5) * (g(1) * (rSges(5,1) * t5 + rSges(5,2) * t4 + t14) + g(2) * (rSges(5,1) * t3 + rSges(5,2) * t2 + t13) + g(3) * (rSges(5,1) * t6 - rSges(5,2) * t25 + t31) + t24) - m(6) * (g(1) * (-t43 * t4 + t14) + g(2) * (-t43 * t2 + t13) + g(3) * (t25 * t43 + t31) + t24 + (g(1) * t5 + g(2) * t3 + g(3) * t6) * (t22 * rSges(6,1) - t20 * rSges(6,2) + pkin(4))), t32 * (-g(3) * t23 + t41), t40 * (-g(1) * t18 + g(2) * t19), -m(6) * (g(1) * ((-t18 * t22 - t20 * t4) * rSges(6,1) + (t18 * t20 - t22 * t4) * rSges(6,2)) + g(2) * ((t19 * t22 - t2 * t20) * rSges(6,1) + (-t19 * t20 - t2 * t22) * rSges(6,2)) - g(3) * (-t20 * rSges(6,1) - t22 * rSges(6,2)) * t25)];
taug = t1(:);
