% Calculate joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:40
% DurationCPUTime: 0.44s
% Computational Cost: add. (1114->92), mult. (674->132), div. (0->0), fcn. (550->10), ass. (0->56)
t68 = qJ(1) + pkin(8);
t64 = qJ(3) + t68;
t57 = sin(t64);
t53 = t57 ^ 2;
t58 = cos(t64);
t54 = t58 ^ 2;
t67 = pkin(9) + qJ(5);
t60 = sin(t67);
t101 = Icges(6,5) * t60;
t62 = cos(t67);
t100 = Icges(6,6) * t62;
t29 = t100 + t101;
t99 = -rSges(6,1) * t62 + rSges(6,2) * t60;
t98 = rSges(5,3) + qJ(4);
t69 = sin(pkin(9));
t70 = cos(pkin(9));
t97 = rSges(5,1) * t70 - rSges(5,2) * t69 + pkin(3);
t96 = t57 * t58;
t35 = t60 * rSges(6,1) + t62 * rSges(6,2);
t93 = m(6) * t35;
t26 = t57 * rSges(4,1) + t58 * rSges(4,2);
t88 = t54 + t53;
t61 = sin(t68);
t72 = sin(qJ(1));
t65 = t72 * pkin(1);
t87 = pkin(2) * t61 + t65;
t63 = cos(t68);
t73 = cos(qJ(1));
t66 = t73 * pkin(1);
t86 = pkin(2) * t63 + t66;
t83 = t29 * t54 + (t101 / 0.2e1 + t100 / 0.2e1 + t29 / 0.2e1) * t53;
t27 = t58 * rSges(4,1) - t57 * rSges(4,2);
t82 = Icges(5,2) * t70 ^ 2 + Icges(6,2) * t62 ^ 2 + Icges(4,3) + (Icges(5,1) * t69 + 0.2e1 * Icges(5,4) * t70) * t69 + (Icges(6,1) * t60 + 0.2e1 * Icges(6,4) * t62) * t60;
t76 = Icges(6,5) * t62 - Icges(6,6) * t60;
t75 = -t57 * rSges(6,3) + t99 * t58;
t74 = -t58 * rSges(6,3) - t99 * t57;
t13 = t98 * t57 + t97 * t58;
t59 = t70 * pkin(4) + pkin(3);
t71 = -pkin(7) - qJ(4);
t10 = t57 * t59 + t58 * t71 + t74;
t11 = -t57 * t71 + t58 * t59 - t75;
t12 = t97 * t57 - t98 * t58;
t44 = t73 * rSges(2,1) - t72 * rSges(2,2);
t43 = t72 * rSges(2,1) + t73 * rSges(2,2);
t25 = t63 * rSges(3,1) - t61 * rSges(3,2) + t66;
t24 = t61 * rSges(3,1) + t63 * rSges(3,2) + t65;
t21 = t27 + t86;
t20 = t87 + t26;
t15 = -Icges(6,3) * t57 - t76 * t58;
t14 = -Icges(6,3) * t58 + t76 * t57;
t9 = t13 + t86;
t8 = t12 + t87;
t7 = t11 + t86;
t6 = t10 + t87;
t3 = t57 * t74 - t58 * t75;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t82; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t26 * t20 + t27 * t21) + m(5) * (t12 * t8 + t13 * t9) + m(6) * (t10 * t6 + t11 * t7) + t82; 0; m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + t82; m(5) * (-t57 * t8 - t58 * t9) + m(6) * (-t57 * t6 - t58 * t7); 0; m(6) * (-t57 * t10 - t58 * t11) + m(5) * (-t57 * t12 - t58 * t13); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t88; (-t57 * t7 + t58 * t6) * t93 + t83; m(6) * t3; (t10 * t58 - t11 * t57) * t93 + t83; 0; m(6) * (t88 * t35 ^ 2 + t3 ^ 2) - t58 * (t54 * t14 + t15 * t96) - t57 * (t14 * t96 + t53 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
