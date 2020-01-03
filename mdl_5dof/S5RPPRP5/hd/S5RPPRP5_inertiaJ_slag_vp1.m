% Calculate joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:23
% DurationCPUTime: 0.62s
% Computational Cost: add. (771->145), mult. (1860->215), div. (0->0), fcn. (1978->6), ass. (0->68)
t107 = Icges(5,1) + Icges(6,1);
t106 = Icges(5,4) - Icges(6,5);
t102 = Icges(6,4) + Icges(5,5);
t105 = Icges(5,2) + Icges(6,3);
t100 = Icges(5,6) - Icges(6,6);
t104 = rSges(6,1) + pkin(4);
t103 = rSges(6,3) + qJ(5);
t101 = Icges(6,2) + Icges(5,3);
t62 = cos(pkin(7));
t63 = sin(qJ(4));
t61 = sin(pkin(7));
t88 = cos(qJ(4));
t75 = t61 * t88;
t46 = -t62 * t63 + t75;
t64 = sin(qJ(1));
t40 = t46 * t64;
t45 = t61 * t63 + t62 * t88;
t41 = t45 * t64;
t65 = cos(qJ(1));
t89 = t100 * t65 + t105 * t40 + t106 * t41;
t86 = t62 * t65;
t42 = t63 * t86 - t65 * t75;
t43 = t45 * t65;
t83 = -t100 * t64 - t105 * t42 + t106 * t43;
t82 = -t102 * t65 - t106 * t40 - t107 * t41;
t81 = t102 * t64 + t106 * t42 - t107 * t43;
t99 = -t100 * t45 + t102 * t46;
t98 = t105 * t45 - t106 * t46;
t97 = -t106 * t45 + t107 * t46;
t96 = t62 ^ 2;
t60 = t65 ^ 2;
t95 = m(4) / 0.2e1;
t94 = m(6) / 0.2e1;
t91 = -rSges(6,2) - pkin(6);
t90 = -rSges(5,3) - pkin(6);
t85 = t100 * t40 + t101 * t65 + t102 * t41;
t84 = t100 * t42 + t101 * t64 - t102 * t43;
t80 = t103 * t45 + t104 * t46;
t79 = t43 * rSges(5,1) - t42 * rSges(5,2);
t78 = t65 * pkin(1) + t64 * qJ(2);
t77 = t64 ^ 2 + t60;
t76 = qJ(3) * t61;
t74 = pkin(2) * t86 + t65 * t76 + t78;
t73 = t95 + m(5) / 0.2e1 + t94;
t72 = pkin(3) * t86 + t74;
t71 = rSges(3,1) * t62 - rSges(3,2) * t61;
t70 = t41 * rSges(5,1) + t40 * rSges(5,2);
t69 = t103 * t42 + t104 * t43;
t56 = t65 * qJ(2);
t66 = t56 + (-t76 - pkin(1) + (-pkin(2) - pkin(3)) * t62) * t64;
t5 = t65 * t90 + t66 - t70;
t6 = t64 * t90 + t72 + t79;
t68 = m(5) * (t5 * t65 + t6 * t64);
t67 = t103 * t40 - t104 * t41;
t48 = rSges(2,1) * t65 - rSges(2,2) * t64;
t47 = -rSges(2,1) * t64 - rSges(2,2) * t65;
t34 = t64 * rSges(3,3) + t65 * t71 + t78;
t33 = t65 * rSges(3,3) + t56 + (-pkin(1) - t71) * t64;
t32 = rSges(5,1) * t46 - rSges(5,2) * t45;
t22 = t64 * rSges(4,2) + (rSges(4,1) * t62 + rSges(4,3) * t61) * t65 + t74;
t21 = t65 * rSges(4,2) + t56 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t62 + (-rSges(4,3) - qJ(3)) * t61) * t64;
t8 = t80 * t65;
t7 = t80 * t64;
t4 = -t64 * t70 - t65 * t79;
t3 = t64 * t91 + t69 + t72;
t2 = t65 * t91 + t66 + t67;
t1 = t64 * t67 - t65 * t69;
t9 = [Icges(2,3) + t97 * t46 + t98 * t45 + (Icges(4,3) + Icges(3,2)) * t96 + ((Icges(4,1) + Icges(3,1)) * t61 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t62) * t61 + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(3) * (t33 ^ 2 + t34 ^ 2) + m(4) * (t21 ^ 2 + t22 ^ 2) + m(2) * (t47 ^ 2 + t48 ^ 2); m(6) * (t2 * t64 - t3 * t65) + m(5) * (t5 * t64 - t6 * t65) + m(3) * (t33 * t64 - t34 * t65) + m(4) * (t21 * t64 - t22 * t65); 0.2e1 * (m(3) / 0.2e1 + t73) * t77; 0.2e1 * ((t2 * t65 + t3 * t64) * t94 + t68 / 0.2e1 + (t21 * t65 + t22 * t64) * t95) * t61; 0; 0.2e1 * t73 * (t61 ^ 2 * t77 + t96); m(6) * (t2 * t8 + t3 * t7) + t32 * t68 - (t98 * t42 + t97 * t43 - t83 * t45 - t81 * t46 - t99 * t64) * t64 / 0.2e1 + (-t98 * t40 + t97 * t41 - t89 * t45 - t82 * t46 + t99 * t65) * t65 / 0.2e1; m(6) * (t64 * t8 - t65 * t7); m(5) * (t32 * t61 * t77 - t4 * t62) + m(6) * (-t1 * t62 + (t64 * t7 + t65 * t8) * t61); m(5) * (t32 ^ 2 * t77 + t4 ^ 2) + m(6) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + (t40 * t89 - t41 * t82 + t65 * t85) * t60 + ((-t42 * t83 - t43 * t81 + t64 * t84) * t64 + (-t40 * t83 + t41 * t81 + t42 * t89 + t43 * t82 + t64 * t85 + t65 * t84) * t65) * t64; m(6) * (t2 * t42 - t3 * t40); m(6) * (t40 * t65 + t42 * t64); m(6) * (-t45 * t62 + (-t40 * t64 + t42 * t65) * t61); m(6) * (t1 * t45 - t40 * t7 + t42 * t8); m(6) * (t40 ^ 2 + t42 ^ 2 + t45 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
