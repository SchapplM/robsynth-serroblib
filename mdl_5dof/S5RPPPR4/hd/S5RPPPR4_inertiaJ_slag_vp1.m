% Calculate joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:07
% DurationCPUTime: 0.34s
% Computational Cost: add. (573->79), mult. (482->108), div. (0->0), fcn. (382->8), ass. (0->43)
t40 = qJ(1) + pkin(7);
t35 = sin(t40);
t32 = t35 ^ 2;
t37 = cos(t40);
t33 = t37 ^ 2;
t23 = t32 + t33;
t39 = pkin(8) + qJ(5);
t34 = sin(t39);
t36 = cos(t39);
t66 = rSges(6,1) * t34 + rSges(6,2) * t36;
t65 = t35 * t37;
t16 = m(6) * t23;
t41 = sin(pkin(8));
t64 = pkin(4) * t41;
t44 = sin(qJ(1));
t63 = t44 * pkin(1);
t60 = m(5) * t23 + t16;
t57 = rSges(5,3) + qJ(4);
t56 = t37 * rSges(6,3) + t66 * t35;
t45 = cos(qJ(1));
t38 = t45 * pkin(1);
t55 = t37 * pkin(2) + t35 * qJ(3) + t38;
t54 = t37 * qJ(3) - t63;
t42 = cos(pkin(8));
t53 = rSges(5,1) * t41 + rSges(5,2) * t42;
t47 = Icges(6,5) * t34 + Icges(6,6) * t36;
t43 = -pkin(6) - qJ(4);
t2 = (t66 + t64) * t37 + (-rSges(6,3) - pkin(2) + t43) * t35 + t54;
t3 = t35 * t64 - t37 * t43 + t55 + t56;
t46 = m(6) * (t35 * t2 - t37 * t3);
t27 = t45 * rSges(2,1) - t44 * rSges(2,2);
t26 = -t44 * rSges(2,1) - t45 * rSges(2,2);
t22 = t36 * rSges(6,1) - t34 * rSges(6,2);
t15 = t37 * rSges(3,1) - t35 * rSges(3,2) + t38;
t14 = -t35 * rSges(3,1) - t37 * rSges(3,2) - t63;
t9 = Icges(6,3) * t35 - t47 * t37;
t8 = Icges(6,3) * t37 + t47 * t35;
t7 = -t37 * rSges(4,2) + t35 * rSges(4,3) + t55;
t6 = t37 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t35 + t54;
t5 = t53 * t35 + t57 * t37 + t55;
t4 = t53 * t37 + (-pkin(2) - t57) * t35 + t54;
t1 = -t35 * t56 + (t35 * rSges(6,3) - t37 * t66) * t37;
t10 = [Icges(5,1) * t42 ^ 2 + Icges(6,1) * t36 ^ 2 + Icges(4,1) + Icges(2,3) + Icges(3,3) + (-0.2e1 * Icges(5,4) * t42 + Icges(5,2) * t41) * t41 + m(2) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t14 ^ 2 + t15 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Icges(6,4) * t36 + Icges(6,2) * t34) * t34; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t35 * t6 - t37 * t7) + m(5) * (t35 * t4 - t37 * t5) + t46; 0; m(4) * t23 + t60; m(5) * (t35 * t5 + t37 * t4) + m(6) * (t37 * t2 + t35 * t3); 0; 0; t60; t22 * t46 + t23 * (Icges(6,5) * t36 - Icges(6,6) * t34); m(6) * t1; t22 * t16; 0; m(6) * (t23 * t22 ^ 2 + t1 ^ 2) + t37 * (t33 * t8 + t9 * t65) + t35 * (t32 * t9 + t8 * t65);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
