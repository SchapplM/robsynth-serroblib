% Calculate joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:56
% DurationCPUTime: 0.51s
% Computational Cost: add. (1114->109), mult. (1276->153), div. (0->0), fcn. (1378->8), ass. (0->61)
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t79 = qJ(3) + pkin(8);
t74 = sin(t79);
t75 = cos(t79);
t31 = -t59 * t75 + t62 * t74;
t94 = t31 ^ 2;
t30 = -t59 * t74 - t62 * t75;
t95 = t30 ^ 2;
t105 = t94 + t95;
t100 = t31 * t30;
t54 = t62 * qJ(2);
t99 = -t59 * pkin(1) + t54;
t57 = sin(qJ(5));
t89 = rSges(6,2) * t57;
t76 = -pkin(4) + t89;
t61 = cos(qJ(3));
t52 = t61 * pkin(3) + pkin(2);
t58 = sin(qJ(3));
t88 = t59 * t58;
t84 = pkin(3) * t88 + t62 * t52;
t60 = cos(qJ(5));
t90 = rSges(6,1) * t60;
t85 = t31 * rSges(6,3) - t30 * t90;
t98 = -t31 * pkin(7) - t76 * t30 - t84 - t85;
t87 = t62 * t58;
t83 = pkin(3) * t87 - t59 * t52;
t86 = t30 * rSges(6,3) + t31 * t90;
t97 = -t30 * pkin(7) + t76 * t31 - t83 - t86;
t96 = t105 * (-Icges(6,5) * t57 - Icges(6,6) * t60);
t39 = -t57 * rSges(6,1) - t60 * rSges(6,2);
t93 = m(6) * t39;
t91 = t59 * pkin(2);
t82 = t62 * pkin(1) + t59 * qJ(2);
t33 = -t62 * t61 - t88;
t34 = t59 * t61 - t87;
t22 = t34 * rSges(4,1) + t33 * rSges(4,2);
t23 = t33 * rSges(4,1) - t34 * rSges(4,2);
t68 = -Icges(6,5) * t60 + Icges(6,6) * t57;
t65 = -t31 * rSges(5,1) + t30 * rSges(5,2) - t83;
t64 = t30 * rSges(5,1) + t31 * rSges(5,2) - t84;
t63 = -t60 ^ 2 * Icges(6,2) - Icges(4,3) - Icges(5,3) + (-Icges(6,1) * t57 - 0.2e1 * Icges(6,4) * t60) * t57;
t55 = t62 * pkin(2);
t41 = t62 * rSges(2,1) - t59 * rSges(2,2);
t40 = -t59 * rSges(2,1) - t62 * rSges(2,2);
t29 = t62 * rSges(3,1) + t59 * rSges(3,3) + t82;
t28 = t62 * rSges(3,3) + t54 + (-rSges(3,1) - pkin(1)) * t59;
t21 = -t23 + t55 + t82;
t20 = t54 + (-pkin(1) - pkin(2)) * t59 - t22;
t19 = t55 + t64;
t18 = t65 - t91;
t13 = Icges(6,3) * t31 + t68 * t30;
t12 = -Icges(6,3) * t30 + t68 * t31;
t11 = -t64 + t82;
t10 = -t65 + t99;
t7 = t55 + t98;
t6 = -t91 + t97;
t5 = t82 - t98;
t4 = -t97 + t99;
t1 = t31 * (t31 * t89 - t86) + t30 * (t30 * t89 + t85);
t2 = [Icges(3,2) + Icges(2,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(2) * (t40 ^ 2 + t41 ^ 2) - t63; m(6) * (t59 * t4 - t62 * t5) + m(5) * (t59 * t10 - t62 * t11) + m(4) * (t59 * t20 - t62 * t21) + m(3) * (t59 * t28 - t62 * t29); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t59 ^ 2 + t62 ^ 2); m(6) * (t6 * t4 + t7 * t5) + m(5) * (t18 * t10 + t19 * t11) + m(4) * (t22 * t20 + t23 * t21) + t63; m(4) * (t22 * t59 - t23 * t62) + m(5) * (t18 * t59 - t19 * t62) + m(6) * (t6 * t59 - t7 * t62); m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) - t63; 0; 0; 0; m(5) + m(6); (-t30 * t4 - t31 * t5) * t93 + t96; (-t30 * t59 + t31 * t62) * t93; (-t30 * t6 - t31 * t7) * t93 - t96; -m(6) * t1; m(6) * (t105 * t39 ^ 2 + t1 ^ 2) + t31 * (-t12 * t100 + t94 * t13) - t30 * (-t13 * t100 + t95 * t12);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
