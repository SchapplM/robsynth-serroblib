% Calculate joint inertia matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:43
% DurationCPUTime: 0.36s
% Computational Cost: add. (934->88), mult. (632->121), div. (0->0), fcn. (502->8), ass. (0->52)
t57 = qJ(1) + pkin(8);
t55 = qJ(3) + t57;
t51 = sin(t55);
t48 = t51 ^ 2;
t52 = cos(t55);
t49 = t52 ^ 2;
t60 = cos(qJ(5));
t88 = Icges(6,5) * t60;
t58 = sin(qJ(5));
t87 = Icges(6,6) * t58;
t33 = -t87 + t88;
t86 = rSges(6,1) * t58 + rSges(6,2) * t60;
t85 = t51 * t52;
t59 = sin(qJ(1));
t82 = t59 * pkin(1);
t79 = t86 * t52;
t78 = t52 * pkin(3) + t51 * qJ(4);
t77 = t48 + t49;
t54 = cos(t57);
t61 = cos(qJ(1));
t56 = t61 * pkin(1);
t76 = pkin(2) * t54 + t56;
t73 = t52 * rSges(6,3) + t51 * t86;
t72 = t33 * t49 + (t88 / 0.2e1 - t87 / 0.2e1 + t33 / 0.2e1) * t48;
t25 = t52 * rSges(4,1) - rSges(4,2) * t51;
t53 = sin(t57);
t71 = -pkin(2) * t53 - t82;
t24 = -rSges(4,1) * t51 - rSges(4,2) * t52;
t40 = t52 * qJ(4);
t12 = t52 * rSges(5,3) + t40 + (rSges(5,2) - pkin(3)) * t51;
t65 = Icges(6,5) * t58 + Icges(6,6) * t60;
t13 = -rSges(5,2) * t52 + t51 * rSges(5,3) + t78;
t9 = t52 * pkin(7) + t73 + t78;
t64 = Icges(6,1) * t60 ^ 2 + Icges(5,1) + Icges(4,3) + (-0.2e1 * Icges(6,4) * t60 + Icges(6,2) * t58) * t58;
t8 = t40 + (-rSges(6,3) - pkin(3) - pkin(7)) * t51 + t79;
t6 = t71 + t8;
t7 = t9 + t76;
t63 = m(6) * (t51 * t6 - t52 * t7);
t62 = m(6) * (t51 * t8 - t52 * t9);
t38 = rSges(2,1) * t61 - rSges(2,2) * t59;
t37 = rSges(6,1) * t60 - rSges(6,2) * t58;
t36 = -rSges(2,1) * t59 - rSges(2,2) * t61;
t23 = rSges(3,1) * t54 - rSges(3,2) * t53 + t56;
t22 = -rSges(3,1) * t53 - rSges(3,2) * t54 - t82;
t21 = t25 + t76;
t20 = t24 + t71;
t15 = Icges(6,3) * t51 - t52 * t65;
t14 = Icges(6,3) * t52 + t51 * t65;
t11 = t13 + t76;
t10 = t12 + t71;
t3 = t52 * (rSges(6,3) * t51 - t79) - t51 * t73;
t1 = [Icges(2,3) + Icges(3,3) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t22 ^ 2 + t23 ^ 2) + m(2) * (t36 ^ 2 + t38 ^ 2) + t64; 0; m(3) + m(4) + m(5) + m(6); m(5) * (t10 * t12 + t11 * t13) + m(6) * (t6 * t8 + t7 * t9) + m(4) * (t20 * t24 + t21 * t25) + t64; 0; m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t64; m(5) * (t10 * t51 - t11 * t52) + t63; 0; m(5) * (t12 * t51 - t13 * t52) + t62; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t77; t37 * t63 + t72; m(6) * t3; t37 * t62 + t72; m(6) * t77 * t37; m(6) * (t37 ^ 2 * t77 + t3 ^ 2) + t52 * (t49 * t14 + t15 * t85) + t51 * (t14 * t85 + t48 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
