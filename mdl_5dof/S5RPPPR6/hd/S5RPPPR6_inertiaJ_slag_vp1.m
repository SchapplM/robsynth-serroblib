% Calculate joint inertia matrix for
% S5RPPPR6
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:39
% DurationCPUTime: 0.50s
% Computational Cost: add. (923->160), mult. (2264->240), div. (0->0), fcn. (2537->8), ass. (0->80)
t65 = sin(pkin(7));
t99 = 0.2e1 * t65;
t60 = t65 ^ 2;
t67 = cos(pkin(7));
t98 = 0.2e1 * t67;
t97 = m(4) / 0.2e1;
t96 = m(5) / 0.2e1;
t95 = m(6) / 0.2e1;
t66 = cos(pkin(8));
t94 = t67 * t66;
t68 = sin(qJ(5));
t93 = t67 * t68;
t70 = cos(qJ(5));
t92 = t67 * t70;
t71 = cos(qJ(1));
t91 = t67 * t71;
t64 = sin(pkin(8));
t69 = sin(qJ(1));
t90 = t69 * t64;
t89 = t69 * t66;
t88 = t71 * t64;
t87 = t71 * t66;
t86 = -pkin(2) - qJ(4);
t85 = t71 * pkin(1) + t69 * qJ(2);
t56 = t71 * qJ(2);
t84 = t71 * pkin(3) + t56;
t83 = t69 ^ 2 + t71 ^ 2;
t82 = qJ(3) * t65;
t81 = t96 + t95;
t42 = t64 * t93 + t65 * t70;
t43 = -t64 * t92 + t65 * t68;
t24 = Icges(6,5) * t43 + Icges(6,6) * t42 + Icges(6,3) * t94;
t25 = Icges(6,4) * t43 + Icges(6,2) * t42 + Icges(6,6) * t94;
t26 = Icges(6,1) * t43 + Icges(6,4) * t42 + Icges(6,5) * t94;
t80 = t24 * t94 + t42 * t25 + t43 * t26;
t45 = t65 * t88 + t89;
t32 = -t45 * t68 + t70 * t91;
t33 = t45 * t70 + t68 * t91;
t44 = -t65 * t87 + t90;
t17 = t33 * rSges(6,1) + t32 * rSges(6,2) + t44 * rSges(6,3);
t79 = -pkin(1) - t82;
t78 = pkin(2) * t91 + t71 * t82 + t85;
t77 = t97 + t81;
t46 = t65 * t89 + t88;
t47 = t65 * t90 - t87;
t34 = -t47 * t68 + t69 * t92;
t35 = t47 * t70 + t69 * t93;
t74 = -t35 * rSges(6,1) - t34 * rSges(6,2);
t18 = -t46 * rSges(6,3) - t74;
t27 = t43 * rSges(6,1) + t42 * rSges(6,2) + rSges(6,3) * t94;
t10 = -t18 * t94 - t46 * t27;
t9 = t17 * t94 - t44 * t27;
t76 = t10 * t71 + t69 * t9;
t75 = rSges(3,1) * t67 - rSges(3,2) * t65;
t73 = t69 * pkin(3) + qJ(4) * t91 + t78;
t19 = -t47 * rSges(5,1) - t46 * rSges(5,2) + ((-rSges(5,3) + t86) * t67 + t79) * t69 + t84;
t20 = t45 * rSges(5,1) - t44 * rSges(5,2) + rSges(5,3) * t91 + t73;
t7 = -t47 * pkin(4) + (rSges(6,3) + pkin(6)) * t46 + (t86 * t67 + t79) * t69 + t74 + t84;
t8 = t45 * pkin(4) + t44 * pkin(6) + t17 + t73;
t72 = (t69 * t8 + t7 * t71) * t95 + (t19 * t71 + t20 * t69) * t96;
t61 = t67 ^ 2;
t49 = t71 * rSges(2,1) - t69 * rSges(2,2);
t48 = -t69 * rSges(2,1) - t71 * rSges(2,2);
t37 = t69 * rSges(3,3) + t75 * t71 + t85;
t36 = t71 * rSges(3,3) + t56 + (-pkin(1) - t75) * t69;
t29 = t69 * rSges(4,1) + (-rSges(4,2) * t67 + rSges(4,3) * t65) * t71 + t78;
t28 = t71 * rSges(4,1) + t56 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t67 + (-rSges(4,3) - qJ(3)) * t65) * t69;
t16 = Icges(6,1) * t35 + Icges(6,4) * t34 - Icges(6,5) * t46;
t15 = Icges(6,1) * t33 + Icges(6,4) * t32 + Icges(6,5) * t44;
t14 = Icges(6,4) * t35 + Icges(6,2) * t34 - Icges(6,6) * t46;
t13 = Icges(6,4) * t33 + Icges(6,2) * t32 + Icges(6,6) * t44;
t12 = Icges(6,5) * t35 + Icges(6,6) * t34 - Icges(6,3) * t46;
t11 = Icges(6,5) * t33 + Icges(6,6) * t32 + Icges(6,3) * t44;
t6 = t80 * t94;
t5 = t46 * t17 + t44 * t18;
t4 = -t46 * t24 + t34 * t25 + t35 * t26;
t3 = t44 * t24 + t32 * t25 + t33 * t26;
t2 = t12 * t94 + t42 * t14 + t43 * t16;
t1 = t11 * t94 + t42 * t13 + t43 * t15;
t21 = [Icges(2,3) + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t60 + m(6) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + m(2) * (t48 ^ 2 + t49 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + ((t66 ^ 2 * Icges(5,2) + Icges(3,2) + Icges(4,3) + (Icges(5,1) * t64 + 0.2e1 * Icges(5,4) * t66) * t64) * t67 + (-Icges(5,5) * t64 - Icges(5,6) * t66 + Icges(3,4) + Icges(4,6)) * t99) * t67 + t80; m(6) * (t69 * t7 - t71 * t8) + m(4) * (t69 * t28 - t71 * t29) + m(5) * (t69 * t19 - t71 * t20) + m(3) * (t69 * t36 - t71 * t37); 0.2e1 * (m(3) / 0.2e1 + t77) * t83; ((t28 * t71 + t29 * t69) * t97 + t72) * t99; 0; 0.2e1 * t77 * (t83 * t60 + t61); t72 * t98; 0; t81 * (-0.1e1 + t83) * t65 * t98; 0.2e1 * t81 * (t83 * t61 + t60); t6 + m(6) * (t10 * t7 + t9 * t8) + (-t2 / 0.2e1 - t4 / 0.2e1) * t46 + (t1 / 0.2e1 + t3 / 0.2e1) * t44; m(6) * (t10 * t69 - t9 * t71); m(6) * (-t5 * t67 + t76 * t65); m(6) * (t5 * t65 + t76 * t67); m(6) * (t10 ^ 2 + t5 ^ 2 + t9 ^ 2) + t44 * ((t44 * t11 + t32 * t13 + t33 * t15) * t44 - (t44 * t12 + t32 * t14 + t33 * t16) * t46 + t3 * t94) - t46 * ((-t46 * t11 + t34 * t13 + t35 * t15) * t44 - (-t46 * t12 + t34 * t14 + t35 * t16) * t46 + t4 * t94) + (t1 * t44 - t2 * t46 + t6) * t94;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
