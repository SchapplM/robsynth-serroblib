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
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:35
% EndTime: 2022-01-23 09:20:35
% DurationCPUTime: 0.40s
% Computational Cost: add. (1697->142), mult. (1430->212), div. (0->0), fcn. (1440->10), ass. (0->71)
t73 = sin(qJ(1));
t93 = t73 * pkin(1);
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t72 = sin(qJ(5));
t74 = cos(qJ(5));
t38 = -Icges(6,3) * t71 + (Icges(6,5) * t74 - Icges(6,6) * t72) * t70;
t39 = -Icges(6,6) * t71 + (Icges(6,4) * t74 - Icges(6,2) * t72) * t70;
t40 = -Icges(6,5) * t71 + (Icges(6,1) * t74 - Icges(6,4) * t72) * t70;
t16 = -t71 * t38 + (-t39 * t72 + t40 * t74) * t70;
t92 = t16 * t71;
t69 = qJ(1) + pkin(8);
t67 = qJ(3) + t69;
t63 = sin(t67);
t91 = t63 * t70;
t64 = cos(t67);
t90 = t64 * t70;
t89 = t64 * t71;
t88 = t71 * t72;
t87 = t71 * t74;
t86 = t64 * pkin(3) + t63 * qJ(4);
t66 = cos(t69);
t75 = cos(qJ(1));
t68 = t75 * pkin(1);
t85 = pkin(2) * t66 + t68;
t84 = Icges(6,5) * t70;
t83 = Icges(6,6) * t70;
t82 = Icges(6,3) * t70;
t36 = t63 * t74 - t64 * t88;
t37 = t63 * t72 + t64 * t87;
t24 = t37 * rSges(6,1) + t36 * rSges(6,2) + rSges(6,3) * t90;
t45 = t64 * rSges(4,1) - t63 * rSges(4,2);
t65 = sin(t69);
t79 = -pkin(2) * t65 - t93;
t44 = -t63 * rSges(4,1) - t64 * rSges(4,2);
t34 = -t63 * t88 - t64 * t74;
t35 = t63 * t87 - t64 * t72;
t78 = -t35 * rSges(6,1) - t34 * rSges(6,2);
t17 = Icges(6,5) * t35 + Icges(6,6) * t34 + t63 * t82;
t19 = Icges(6,4) * t35 + Icges(6,2) * t34 + t63 * t83;
t21 = Icges(6,1) * t35 + Icges(6,4) * t34 + t63 * t84;
t3 = -t71 * t17 + (-t19 * t72 + t21 * t74) * t70;
t18 = Icges(6,5) * t37 + Icges(6,6) * t36 + t64 * t82;
t20 = Icges(6,4) * t37 + Icges(6,2) * t36 + t64 * t83;
t22 = Icges(6,1) * t37 + Icges(6,4) * t36 + t64 * t84;
t4 = -t71 * t18 + (-t20 * t72 + t22 * t74) * t70;
t8 = t34 * t39 + t35 * t40 + t38 * t91;
t9 = t36 * t39 + t37 * t40 + t38 * t90;
t77 = -t92 + (t3 + t8) * t91 / 0.2e1 + (t4 + t9) * t90 / 0.2e1;
t13 = pkin(4) * t89 + pkin(7) * t90 + t24 + t86;
t28 = rSges(5,1) * t89 - rSges(5,2) * t90 + t63 * rSges(5,3) + t86;
t57 = t64 * qJ(4);
t27 = rSges(5,2) * t91 + t57 + t64 * rSges(5,3) + (-rSges(5,1) * t71 - pkin(3)) * t63;
t76 = Icges(5,2) * t71 ^ 2 + Icges(4,3) + t16 + (Icges(5,1) * t70 + 0.2e1 * Icges(5,4) * t71) * t70;
t12 = t57 + (-pkin(4) * t71 - pkin(3) + (-rSges(6,3) - pkin(7)) * t70) * t63 + t78;
t55 = t75 * rSges(2,1) - t73 * rSges(2,2);
t54 = -t73 * rSges(2,1) - t75 * rSges(2,2);
t43 = t66 * rSges(3,1) - t65 * rSges(3,2) + t68;
t42 = -t65 * rSges(3,1) - t66 * rSges(3,2) - t93;
t41 = -t71 * rSges(6,3) + (rSges(6,1) * t74 - rSges(6,2) * t72) * t70;
t30 = t45 + t85;
t29 = t44 + t79;
t26 = t28 + t85;
t25 = t27 + t79;
t23 = rSges(6,3) * t91 - t78;
t15 = -t71 * t24 - t41 * t90;
t14 = t71 * t23 + t41 * t91;
t11 = t13 + t85;
t10 = t12 + t79;
t5 = (t23 * t64 - t24 * t63) * t70;
t1 = [Icges(2,3) + Icges(3,3) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(2) * (t54 ^ 2 + t55 ^ 2) + t76; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t12 * t10 + t13 * t11) + m(5) * (t27 * t25 + t28 * t26) + m(4) * (t44 * t29 + t45 * t30) + t76; 0; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t76; m(6) * (t63 * t10 - t64 * t11) + m(5) * (t63 * t25 - t64 * t26); 0; m(6) * (t63 * t12 - t64 * t13) + m(5) * (t63 * t27 - t64 * t28); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t63 ^ 2 + t64 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t77; m(6) * t5; m(6) * (t14 * t12 + t15 * t13) + t77; m(6) * (t14 * t63 - t15 * t64); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) + ((t18 * t90 + t36 * t20 + t37 * t22) * t90 + (t17 * t90 + t36 * t19 + t37 * t21) * t91 - t9 * t71) * t90 + ((t18 * t91 + t34 * t20 + t35 * t22) * t90 + (t17 * t91 + t34 * t19 + t35 * t21) * t91 - t8 * t71) * t91 - t71 * (-t92 + (t3 * t63 + t4 * t64) * t70);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
