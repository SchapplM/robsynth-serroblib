% Calculate joint inertia matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:27
% DurationCPUTime: 0.38s
% Computational Cost: add. (759->81), mult. (628->123), div. (0->0), fcn. (514->8), ass. (0->50)
t58 = qJ(1) + qJ(2);
t54 = sin(t58);
t50 = t54 ^ 2;
t55 = cos(t58);
t51 = t55 ^ 2;
t57 = pkin(7) + qJ(4);
t52 = sin(t57);
t89 = Icges(5,5) * t52;
t53 = cos(t57);
t88 = Icges(5,6) * t53;
t25 = t88 + t89;
t87 = qJ(3) + rSges(4,3);
t59 = sin(pkin(7));
t60 = cos(pkin(7));
t86 = -rSges(4,1) * t60 + rSges(4,2) * t59 - pkin(2);
t85 = t54 * t55;
t30 = t52 * rSges(5,1) + t53 * rSges(5,2);
t82 = m(5) * t30;
t62 = sin(qJ(1));
t81 = t62 * pkin(1);
t79 = rSges(5,1) * t53;
t77 = rSges(5,2) * t52;
t76 = t55 * rSges(5,3) + t54 * t77;
t75 = t50 + t51;
t72 = t25 * t51 + (t89 / 0.2e1 + t88 / 0.2e1 + t25 / 0.2e1) * t50;
t32 = t55 * rSges(3,1) - t54 * rSges(3,2);
t71 = Icges(4,2) * t60 ^ 2 + Icges(5,2) * t53 ^ 2 + Icges(3,3) + (Icges(4,1) * t59 + 0.2e1 * Icges(4,4) * t60) * t59 + (Icges(5,1) * t52 + 0.2e1 * Icges(5,4) * t53) * t52;
t31 = -t54 * rSges(3,1) - t55 * rSges(3,2);
t65 = Icges(5,5) * t53 - Icges(5,6) * t52;
t64 = t54 * rSges(5,3) + (-t77 + t79) * t55;
t13 = t54 * t87 - t86 * t55;
t12 = t86 * t54 + t87 * t55;
t48 = t60 * pkin(3) + pkin(2);
t61 = -pkin(6) - qJ(3);
t9 = t55 * t48 - t54 * t61 + t64;
t8 = -t55 * t61 + (-t48 - t79) * t54 + t76;
t63 = cos(qJ(1));
t56 = t63 * pkin(1);
t38 = t63 * rSges(2,1) - t62 * rSges(2,2);
t37 = -t62 * rSges(2,1) - t63 * rSges(2,2);
t23 = t32 + t56;
t22 = t31 - t81;
t15 = Icges(5,3) * t54 + t65 * t55;
t14 = -Icges(5,3) * t55 + t65 * t54;
t11 = t13 + t56;
t10 = t12 - t81;
t7 = t56 + t9;
t6 = t8 - t81;
t5 = t54 * (t54 * t79 - t76) + t55 * t64;
t1 = [Icges(2,3) + m(2) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t71; m(3) * (t31 * t22 + t32 * t23) + m(4) * (t12 * t10 + t13 * t11) + m(5) * (t8 * t6 + t9 * t7) + t71; m(3) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t71; m(4) * (t54 * t10 - t55 * t11) + m(5) * (t54 * t6 - t55 * t7); m(4) * (t54 * t12 - t55 * t13) + m(5) * (t54 * t8 - t55 * t9); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t75; (-t54 * t7 - t55 * t6) * t82 + t72; (-t54 * t9 - t55 * t8) * t82 + t72; 0; m(5) * (t75 * t30 ^ 2 + t5 ^ 2) + t54 * (-t14 * t85 + t50 * t15) - t55 * (t51 * t14 - t15 * t85);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
