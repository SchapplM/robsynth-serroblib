% Calculate joint inertia matrix for
% S5RPRPR3
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:57
% EndTime: 2020-01-03 11:36:00
% DurationCPUTime: 0.54s
% Computational Cost: add. (1697->141), mult. (1430->212), div. (0->0), fcn. (1440->10), ass. (0->71)
t76 = sin(pkin(9));
t77 = cos(pkin(9));
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t38 = -Icges(6,3) * t77 + (Icges(6,5) * t80 - Icges(6,6) * t78) * t76;
t39 = -Icges(6,6) * t77 + (Icges(6,4) * t80 - Icges(6,2) * t78) * t76;
t40 = -Icges(6,5) * t77 + (Icges(6,1) * t80 - Icges(6,4) * t78) * t76;
t16 = -t77 * t38 + (-t39 * t78 + t40 * t80) * t76;
t98 = t16 * t77;
t75 = qJ(1) + pkin(8);
t72 = qJ(3) + t75;
t68 = sin(t72);
t97 = t68 * t76;
t96 = t68 * t77;
t69 = cos(t72);
t95 = t69 * t76;
t94 = t69 * t77;
t93 = t77 * t78;
t92 = t77 * t80;
t91 = t69 * pkin(3) + t68 * qJ(4);
t44 = t68 * rSges(4,1) + t69 * rSges(4,2);
t70 = sin(t75);
t79 = sin(qJ(1));
t73 = t79 * pkin(1);
t90 = pkin(2) * t70 + t73;
t71 = cos(t75);
t81 = cos(qJ(1));
t74 = t81 * pkin(1);
t89 = pkin(2) * t71 + t74;
t88 = Icges(6,5) * t76;
t87 = Icges(6,6) * t76;
t86 = Icges(6,3) * t76;
t34 = -t68 * t93 - t69 * t80;
t35 = t68 * t92 - t69 * t78;
t23 = t35 * rSges(6,1) + t34 * rSges(6,2) + rSges(6,3) * t97;
t45 = t69 * rSges(4,1) - t68 * rSges(4,2);
t17 = Icges(6,5) * t35 + Icges(6,6) * t34 + t68 * t86;
t19 = Icges(6,4) * t35 + Icges(6,2) * t34 + t68 * t87;
t21 = Icges(6,1) * t35 + Icges(6,4) * t34 + t68 * t88;
t3 = -t77 * t17 + (-t19 * t78 + t21 * t80) * t76;
t36 = -t68 * t80 + t69 * t93;
t37 = -t68 * t78 - t69 * t92;
t18 = Icges(6,5) * t37 + Icges(6,6) * t36 - t69 * t86;
t20 = Icges(6,4) * t37 + Icges(6,2) * t36 - t69 * t87;
t22 = Icges(6,1) * t37 + Icges(6,4) * t36 - t69 * t88;
t4 = -t77 * t18 + (-t20 * t78 + t22 * t80) * t76;
t8 = t34 * t39 + t35 * t40 + t38 * t97;
t9 = t36 * t39 + t37 * t40 - t38 * t95;
t83 = -t98 + (t3 + t8) * t97 / 0.2e1 - (t4 + t9) * t95 / 0.2e1;
t24 = t37 * rSges(6,1) + t36 * rSges(6,2) - rSges(6,3) * t95;
t28 = rSges(5,1) * t94 - rSges(5,2) * t95 + t68 * rSges(5,3) + t91;
t64 = t68 * pkin(3);
t12 = pkin(4) * t96 + pkin(7) * t97 - t69 * qJ(4) + t23 + t64;
t27 = -rSges(5,2) * t97 + rSges(5,1) * t96 + t64 + (-rSges(5,3) - qJ(4)) * t69;
t13 = pkin(4) * t94 + pkin(7) * t95 - t24 + t91;
t82 = Icges(5,2) * t77 ^ 2 + Icges(4,3) + t16 + (Icges(5,1) * t76 + 0.2e1 * Icges(5,4) * t77) * t76;
t58 = t81 * rSges(2,1) - t79 * rSges(2,2);
t57 = t79 * rSges(2,1) + t81 * rSges(2,2);
t43 = t71 * rSges(3,1) - t70 * rSges(3,2) + t74;
t42 = t70 * rSges(3,1) + t71 * rSges(3,2) + t73;
t41 = -t77 * rSges(6,3) + (rSges(6,1) * t80 - rSges(6,2) * t78) * t76;
t30 = t45 + t89;
t29 = t90 + t44;
t26 = t28 + t89;
t25 = t27 + t90;
t15 = t77 * t24 - t41 * t95;
t14 = -t77 * t23 - t41 * t97;
t11 = t13 + t89;
t10 = t12 + t90;
t5 = (t23 * t69 + t24 * t68) * t76;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t57 ^ 2 + t58 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t82; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t44 * t29 + t45 * t30) + m(5) * (t27 * t25 + t28 * t26) + m(6) * (t12 * t10 + t13 * t11) + t82; 0; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t82; m(5) * (-t68 * t25 - t69 * t26) + m(6) * (-t68 * t10 - t69 * t11); 0; m(6) * (-t68 * t12 - t69 * t13) + m(5) * (-t68 * t27 - t69 * t28); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t68 ^ 2 + t69 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t83; m(6) * t5; m(6) * (t14 * t12 + t15 * t13) + t83; m(6) * (-t14 * t68 - t15 * t69); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) - t77 * (-t98 + (t3 * t68 - t4 * t69) * t76) + (-t8 * t77 + (t17 * t97 + t34 * t19 + t35 * t21) * t97 - (t18 * t97 + t34 * t20 + t35 * t22) * t95) * t97 - (-t9 * t77 + (-t17 * t95 + t36 * t19 + t37 * t21) * t97 - (-t18 * t95 + t36 * t20 + t37 * t22) * t95) * t95;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
