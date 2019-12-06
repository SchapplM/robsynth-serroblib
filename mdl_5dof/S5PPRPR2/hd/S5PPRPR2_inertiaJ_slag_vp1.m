% Calculate joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:02
% DurationCPUTime: 0.67s
% Computational Cost: add. (1410->125), mult. (1908->218), div. (0->0), fcn. (1990->6), ass. (0->70)
t62 = pkin(8) + qJ(3);
t58 = sin(t62);
t59 = cos(t62);
t101 = (-Icges(5,4) + Icges(4,5)) * t59 + (Icges(5,5) - Icges(4,6)) * t58;
t100 = Icges(5,1) + Icges(4,3);
t63 = sin(pkin(7));
t60 = t63 ^ 2;
t64 = cos(pkin(7));
t61 = t64 ^ 2;
t85 = t60 + t61;
t99 = t100 * t64 - t101 * t63;
t98 = t100 * t63 + t101 * t64;
t97 = t58 ^ 2;
t96 = t63 / 0.2e1;
t95 = -m(5) - m(6);
t65 = sin(qJ(5));
t66 = cos(qJ(5));
t29 = Icges(6,3) * t58 + (-Icges(6,5) * t65 - Icges(6,6) * t66) * t59;
t94 = t58 * t29;
t93 = t59 * t63;
t92 = t59 * t64;
t91 = t63 * t65;
t90 = t63 * t66;
t89 = t64 * t65;
t88 = t64 * t66;
t87 = t85 * (pkin(3) * t59 + qJ(4) * t58);
t53 = t58 * pkin(3) - t59 * qJ(4);
t86 = t58 * rSges(5,2) + t59 * rSges(5,3) - t53;
t84 = Icges(6,5) * t59;
t83 = Icges(6,6) * t59;
t82 = Icges(6,3) * t59;
t81 = m(5) / 0.2e1 + m(6) / 0.2e1;
t32 = t58 * rSges(6,3) + (-rSges(6,1) * t65 - rSges(6,2) * t66) * t59;
t80 = -pkin(6) * t58 - t32 - t53;
t55 = t58 * rSges(4,1) + t59 * rSges(4,2);
t51 = t58 * t91 - t88;
t50 = t58 * t90 + t89;
t49 = t58 * t89 + t90;
t48 = t58 * t88 - t91;
t31 = Icges(6,5) * t58 + (-Icges(6,1) * t65 - Icges(6,4) * t66) * t59;
t30 = Icges(6,6) * t58 + (-Icges(6,4) * t65 - Icges(6,2) * t66) * t59;
t28 = t86 * t64;
t27 = t86 * t63;
t26 = t51 * rSges(6,1) + t50 * rSges(6,2) + rSges(6,3) * t93;
t25 = t49 * rSges(6,1) + t48 * rSges(6,2) + rSges(6,3) * t92;
t24 = Icges(6,1) * t51 + Icges(6,4) * t50 + t63 * t84;
t23 = Icges(6,1) * t49 + Icges(6,4) * t48 + t64 * t84;
t22 = Icges(6,4) * t51 + Icges(6,2) * t50 + t63 * t83;
t21 = Icges(6,4) * t49 + Icges(6,2) * t48 + t64 * t83;
t20 = Icges(6,5) * t51 + Icges(6,6) * t50 + t63 * t82;
t19 = Icges(6,5) * t49 + Icges(6,6) * t48 + t64 * t82;
t18 = t85 * (rSges(4,1) * t59 - rSges(4,2) * t58);
t17 = t80 * t64;
t16 = t80 * t63;
t15 = t58 * t25 - t32 * t92;
t14 = -t58 * t26 + t32 * t93;
t13 = t87 + t85 * (-rSges(5,2) * t59 + rSges(5,3) * t58);
t12 = (-t25 * t63 + t26 * t64) * t59;
t11 = t58 * t20 + (-t22 * t66 - t24 * t65) * t59;
t10 = t58 * t19 + (-t21 * t66 - t23 * t65) * t59;
t9 = t85 * t59 * pkin(6) + t64 * t25 + t63 * t26 + t87;
t8 = t20 * t93 + t50 * t22 + t51 * t24;
t7 = t19 * t93 + t50 * t21 + t51 * t23;
t6 = t20 * t92 + t48 * t22 + t49 * t24;
t5 = t19 * t92 + t48 * t21 + t49 * t23;
t4 = t7 * t63 - t8 * t64;
t3 = t5 * t63 - t6 * t64;
t2 = (t50 * t30 + t51 * t31) * t58 + (t7 * t64 + (t8 + t94) * t63) * t59;
t1 = (t48 * t30 + t49 * t31) * t58 + (t6 * t63 + (t5 + t94) * t64) * t59;
t33 = [m(2) + m(3) + m(4) - t95; 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t81) * t85; m(4) * t18 + m(5) * t13 + m(6) * t9; m(5) * (-t27 * t64 + t28 * t63) + m(6) * (-t16 * t64 + t17 * t63); m(4) * (t55 ^ 2 * t85 + t18 ^ 2) + m(5) * (t13 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2 + t9 ^ 2) + (t99 * t61 - t4) * t64 + (t3 + t98 * t60 + (t99 * t63 + t98 * t64) * t64) * t63; t95 * t59; 0; m(5) * (-t59 * t13 + (t27 * t63 + t28 * t64) * t58) + m(6) * (-t59 * t9 + (t16 * t63 + t17 * t64) * t58); 0.2e1 * t81 * (t59 ^ 2 + t85 * t97); m(6) * t12; m(6) * (t14 * t63 - t15 * t64); t58 * (t10 * t63 - t11 * t64) / 0.2e1 + t1 * t96 - t64 * t2 / 0.2e1 + m(6) * (t12 * t9 + t14 * t17 + t15 * t16) + (t64 * t3 / 0.2e1 + t4 * t96) * t59; m(6) * (-t12 * t59 + (t14 * t64 + t15 * t63) * t58); m(6) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) + t1 * t92 + t2 * t93 + t58 * (t97 * t29 + (t10 * t64 + t11 * t63 + (-t30 * t66 - t31 * t65) * t58) * t59);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t33(1), t33(2), t33(4), t33(7), t33(11); t33(2), t33(3), t33(5), t33(8), t33(12); t33(4), t33(5), t33(6), t33(9), t33(13); t33(7), t33(8), t33(9), t33(10), t33(14); t33(11), t33(12), t33(13), t33(14), t33(15);];
Mq = res;
