% Calculate joint inertia matrix for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:49
% DurationCPUTime: 0.27s
% Computational Cost: add. (322->56), mult. (346->83), div. (0->0), fcn. (276->6), ass. (0->33)
t33 = qJ(1) + pkin(6);
t30 = sin(t33);
t28 = t30 ^ 2;
t31 = cos(t33);
t29 = t31 ^ 2;
t50 = t28 + t29;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t56 = rSges(5,1) * t34 + rSges(5,2) * t36;
t55 = t30 * t31;
t54 = t56 * t31;
t35 = sin(qJ(1));
t53 = t35 * pkin(1);
t47 = t31 * rSges(5,3) + t30 * t56;
t37 = cos(qJ(1));
t32 = t37 * pkin(1);
t46 = t31 * pkin(2) + t30 * qJ(3) + t32;
t45 = t31 * qJ(3) - t53;
t39 = Icges(5,5) * t34 + Icges(5,6) * t36;
t2 = t54 + (-rSges(5,3) - pkin(2) - pkin(5)) * t30 + t45;
t3 = t31 * pkin(5) + t46 + t47;
t38 = m(5) * (t30 * t2 - t31 * t3);
t21 = t37 * rSges(2,1) - t35 * rSges(2,2);
t20 = t36 * rSges(5,1) - t34 * rSges(5,2);
t19 = -t35 * rSges(2,1) - t37 * rSges(2,2);
t13 = t31 * rSges(3,1) - t30 * rSges(3,2) + t32;
t12 = -t30 * rSges(3,1) - t31 * rSges(3,2) - t53;
t7 = Icges(5,3) * t30 - t39 * t31;
t6 = Icges(5,3) * t31 + t39 * t30;
t5 = -t31 * rSges(4,2) + t30 * rSges(4,3) + t46;
t4 = t31 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t30 + t45;
t1 = -t30 * t47 + (t30 * rSges(5,3) - t54) * t31;
t8 = [Icges(5,1) * t36 ^ 2 + Icges(4,1) + Icges(2,3) + Icges(3,3) + m(2) * (t19 ^ 2 + t21 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Icges(5,4) * t36 + Icges(5,2) * t34) * t34; 0; m(3) + m(4) + m(5); m(4) * (t30 * t4 - t31 * t5) + t38; 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t50; t20 * t38 + t50 * (Icges(5,5) * t36 - Icges(5,6) * t34); m(5) * t1; m(5) * t50 * t20; m(5) * (t50 * t20 ^ 2 + t1 ^ 2) + t31 * (t29 * t6 + t7 * t55) + t30 * (t28 * t7 + t6 * t55);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
