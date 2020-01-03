% Calculate joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:00
% DurationCPUTime: 0.67s
% Computational Cost: add. (704->105), mult. (920->148), div. (0->0), fcn. (762->6), ass. (0->57)
t56 = pkin(7) + qJ(4);
t50 = cos(t56);
t114 = t50 ^ 2;
t62 = sin(qJ(1));
t57 = t62 ^ 2;
t63 = cos(qJ(1));
t58 = t63 ^ 2;
t44 = t57 + t58;
t112 = Icges(6,4) + Icges(5,5);
t111 = Icges(5,6) - Icges(6,6);
t49 = sin(t56);
t109 = t111 * t50 + t112 * t49;
t103 = rSges(6,1) + pkin(4);
t86 = rSges(6,3) + qJ(5);
t81 = t86 * t50;
t107 = (t103 * t49 - t81) * t63;
t106 = Icges(6,2) + Icges(5,3);
t102 = (rSges(5,1) * t49 + rSges(5,2) * t50) * t63;
t101 = t106 * t63 + t109 * t62;
t100 = t106 * t62 - t109 * t63;
t36 = m(5) * t44;
t97 = m(6) * t50;
t59 = sin(pkin(7));
t96 = pkin(3) * t59;
t95 = t49 * t62;
t94 = t50 * t62;
t93 = t103 * t50 + t86 * t49;
t92 = t63 * pkin(1) + t62 * qJ(2);
t87 = rSges(4,3) + qJ(3);
t85 = t36 + (m(4) + m(6)) * t44;
t84 = rSges(5,1) * t95 + rSges(5,2) * t94 + t63 * rSges(5,3);
t83 = -t63 * rSges(6,2) - t103 * t95;
t52 = t63 * qJ(2);
t61 = -pkin(6) - qJ(3);
t82 = t62 * t61 + t63 * t96 + t52;
t2 = (-rSges(6,2) - pkin(1)) * t62 + t107 + t82;
t65 = -t63 * t61 + t62 * t96 + t92;
t3 = -t62 * t81 + t65 - t83;
t79 = t62 * t2 - t63 * t3;
t7 = t93 * t62;
t8 = t93 * t63;
t78 = t7 * t62 + t8 * t63;
t60 = cos(pkin(7));
t77 = rSges(4,1) * t59 + rSges(4,2) * t60;
t5 = t102 + (-rSges(5,3) - pkin(1)) * t62 + t82;
t6 = t65 + t84;
t64 = m(5) * (t62 * t5 - t63 * t6);
t39 = t63 * rSges(2,1) - t62 * rSges(2,2);
t38 = -t62 * rSges(2,1) - t63 * rSges(2,2);
t34 = t50 * rSges(5,1) - t49 * rSges(5,2);
t24 = -t63 * rSges(3,2) + t62 * rSges(3,3) + t92;
t23 = t63 * rSges(3,3) + t52 + (rSges(3,2) - pkin(1)) * t62;
t10 = t62 * t77 + t63 * t87 + t92;
t9 = t52 + t77 * t63 + (-pkin(1) - t87) * t62;
t4 = -t62 * t84 + (t62 * rSges(5,3) - t102) * t63;
t1 = (t86 * t94 + t83) * t62 + (t62 * rSges(6,2) - t107) * t63;
t11 = [Icges(4,1) * t60 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t60 + Icges(4,2) * t59) * t59 + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t10 ^ 2 + t9 ^ 2) + m(2) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t23 ^ 2 + t24 ^ 2) + (Icges(5,1) + Icges(6,1)) * t114 + (0.2e1 * (Icges(6,5) - Icges(5,4)) * t50 + (Icges(5,2) + Icges(6,3)) * t49) * t49; m(6) * t79 + t64 + m(4) * (-t63 * t10 + t62 * t9) + m(3) * (t62 * t23 - t63 * t24); m(3) * t44 + t85; m(6) * (t63 * t2 + t62 * t3) + m(5) * (t63 * t5 + t62 * t6) + m(4) * (t62 * t10 + t63 * t9); 0; t85; m(6) * (t7 * t2 - t8 * t3) + t34 * t64 + t44 * (-t111 * t49 + t112 * t50); m(6) * t78 + t34 * t36; m(6) * (-t8 * t62 + t7 * t63); m(5) * (t34 ^ 2 * t44 + t4 ^ 2) + m(6) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t100 * t62 * t57 + (t101 * t58 + (t100 * t63 + t101 * t62) * t62) * t63; -t79 * t97; -t44 * t97; 0; m(6) * (t49 * t1 - t50 * t78); m(6) * (t44 * t114 + t49 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
