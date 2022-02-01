% Calculate joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:51
% EndTime: 2022-01-23 09:27:51
% DurationCPUTime: 0.59s
% Computational Cost: add. (1385->110), mult. (970->147), div. (0->0), fcn. (820->8), ass. (0->65)
t130 = Icges(5,5) + Icges(6,5);
t129 = Icges(5,6) + Icges(6,6);
t79 = sin(qJ(4));
t127 = t130 * t79;
t81 = cos(qJ(4));
t128 = t129 * t81;
t122 = t127 + t128;
t77 = qJ(1) + pkin(8);
t75 = qJ(3) + t77;
t70 = sin(t75);
t67 = t70 ^ 2;
t71 = cos(t75);
t68 = t71 ^ 2;
t126 = -t129 * t79 + t130 * t81;
t123 = Icges(5,3) + Icges(6,3);
t120 = t123 * t71 - t126 * t70;
t119 = t123 * t70 + t126 * t71;
t56 = t79 * rSges(5,1) + t81 * rSges(5,2);
t116 = m(5) * t56;
t80 = sin(qJ(1));
t115 = t80 * pkin(1);
t114 = rSges(5,1) * t81;
t113 = t70 * t79;
t112 = t71 * t79;
t111 = t71 * t81;
t110 = -rSges(6,2) * t113 - t71 * rSges(6,3);
t109 = rSges(5,2) * t113 + t71 * rSges(5,3);
t108 = t71 * pkin(3) + t70 * pkin(7);
t107 = t67 + t68;
t74 = cos(t77);
t82 = cos(qJ(1));
t76 = t82 * pkin(1);
t106 = pkin(2) * t74 + t76;
t101 = -t81 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t79;
t72 = t81 * pkin(4) + pkin(3);
t100 = -rSges(6,1) * t81 - t72;
t38 = t71 * rSges(4,1) - t70 * rSges(4,2);
t73 = sin(t77);
t99 = -pkin(2) * t73 - t115;
t98 = Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t81 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t79 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t81) * t79;
t37 = -t70 * rSges(4,1) - t71 * rSges(4,2);
t85 = rSges(5,1) * t111 - rSges(5,2) * t112 + t70 * rSges(5,3);
t84 = rSges(6,1) * t111 - rSges(6,2) * t112 + t70 * rSges(6,3) + t71 * t72;
t83 = t122 * t68 + (t128 / 0.2e1 + t127 / 0.2e1 + t122 / 0.2e1) * t67;
t18 = t85 + t108;
t65 = t71 * pkin(7);
t17 = t65 + (-pkin(3) - t114) * t70 + t109;
t78 = -qJ(5) - pkin(7);
t16 = -t70 * t78 + t84;
t15 = t100 * t70 - t71 * t78 - t110;
t58 = t82 * rSges(2,1) - t80 * rSges(2,2);
t57 = -t80 * rSges(2,1) - t82 * rSges(2,2);
t36 = t74 * rSges(3,1) - t73 * rSges(3,2) + t76;
t35 = -t73 * rSges(3,1) - t74 * rSges(3,2) - t115;
t34 = t38 + t106;
t33 = t37 + t99;
t32 = t101 * t71;
t31 = t101 * t70;
t14 = t18 + t106;
t13 = t17 + t99;
t12 = t16 + t106;
t11 = t15 + t99;
t6 = t70 * (t70 * t114 - t109) + t71 * t85;
t1 = (t84 - t108) * t71 + (t65 + (-pkin(3) - t100) * t70 + t110) * t70;
t2 = [Icges(2,3) + Icges(3,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(2) * (t57 ^ 2 + t58 ^ 2) + t98; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t15 * t11 + t16 * t12) + m(5) * (t17 * t13 + t18 * t14) + m(4) * (t37 * t33 + t38 * t34) + t98; 0; m(6) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + t98; m(6) * (t32 * t11 + t31 * t12) + (-t13 * t71 - t14 * t70) * t116 + t83; m(5) * t6 + m(6) * t1; m(6) * (t32 * t15 + t31 * t16) + (-t17 * t71 - t18 * t70) * t116 + t83; m(5) * (t107 * t56 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t119 * t70 * t67 + (t120 * t68 + (t119 * t71 + t120 * t70) * t70) * t71; m(6) * (t70 * t11 - t71 * t12); 0; m(6) * (t70 * t15 - t71 * t16); m(6) * (-t71 * t31 + t70 * t32); m(6) * t107;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
