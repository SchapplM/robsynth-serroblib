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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:09
% DurationCPUTime: 0.65s
% Computational Cost: add. (1385->111), mult. (970->147), div. (0->0), fcn. (820->8), ass. (0->64)
t130 = Icges(5,5) + Icges(6,5);
t129 = -Icges(5,6) - Icges(6,6);
t85 = cos(qJ(4));
t125 = t129 * t85;
t83 = sin(qJ(4));
t126 = t130 * t83;
t124 = t125 - t126;
t81 = qJ(1) + pkin(8);
t78 = qJ(3) + t81;
t73 = sin(t78);
t69 = t73 ^ 2;
t74 = cos(t78);
t70 = t74 ^ 2;
t128 = t129 * t83 + t130 * t85;
t127 = Icges(5,3) + Icges(6,3);
t122 = t127 * t74 - t128 * t73;
t121 = t127 * t73 + t128 * t74;
t58 = rSges(5,1) * t83 + rSges(5,2) * t85;
t118 = m(5) * t58;
t117 = t73 * t83;
t116 = t73 * t85;
t115 = t74 * t83;
t114 = t74 * t85;
t37 = rSges(4,1) * t73 + rSges(4,2) * t74;
t113 = pkin(3) * t74 + pkin(7) * t73;
t112 = t70 + t69;
t76 = sin(t81);
t84 = sin(qJ(1));
t79 = t84 * pkin(1);
t111 = pkin(2) * t76 + t79;
t77 = cos(t81);
t86 = cos(qJ(1));
t80 = t86 * pkin(1);
t110 = pkin(2) * t77 + t80;
t105 = t85 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t83;
t38 = rSges(4,1) * t74 - rSges(4,2) * t73;
t104 = rSges(5,1) * t116 - rSges(5,2) * t117;
t103 = Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t85 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t83 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t85) * t83;
t90 = -rSges(5,1) * t114 + rSges(5,2) * t115 - rSges(5,3) * t73;
t75 = pkin(4) * t85 + pkin(3);
t82 = -qJ(5) - pkin(7);
t89 = rSges(6,1) * t116 - rSges(6,2) * t117 + t73 * t75 + t74 * t82;
t88 = rSges(6,1) * t114 - rSges(6,2) * t115 + rSges(6,3) * t73 + t74 * t75;
t87 = -t124 * t70 + (-t125 / 0.2e1 + t126 / 0.2e1 - t124 / 0.2e1) * t69;
t18 = -t90 + t113;
t15 = -rSges(6,3) * t74 + t89;
t16 = -t73 * t82 + t88;
t67 = t73 * pkin(3);
t17 = t67 + (-rSges(5,3) - pkin(7)) * t74 + t104;
t60 = rSges(2,1) * t86 - rSges(2,2) * t84;
t59 = rSges(2,1) * t84 + rSges(2,2) * t86;
t36 = rSges(3,1) * t77 - rSges(3,2) * t76 + t80;
t35 = rSges(3,1) * t76 + rSges(3,2) * t77 + t79;
t34 = t38 + t110;
t33 = t111 + t37;
t32 = t105 * t74;
t31 = t105 * t73;
t14 = t18 + t110;
t13 = t17 + t111;
t12 = t16 + t110;
t11 = t15 + t111;
t6 = -t74 * t90 + t73 * (-rSges(5,3) * t74 + t104);
t1 = (t88 - t113) * t74 + (-t67 + (-rSges(6,3) + pkin(7) - t82) * t74 + t89) * t73;
t2 = [Icges(2,3) + Icges(3,3) + m(2) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2) + t103; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t33 * t37 + t34 * t38) + m(5) * (t13 * t17 + t14 * t18) + m(6) * (t11 * t15 + t12 * t16) + t103; 0; m(6) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + t103; m(6) * (t11 * t32 - t12 * t31) + (t13 * t74 - t14 * t73) * t118 + t87; m(5) * t6 + m(6) * t1; m(6) * (t15 * t32 - t16 * t31) + (t17 * t74 - t18 * t73) * t118 + t87; m(5) * (t112 * t58 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t122 * t74 * t70 + (t121 * t69 + (t121 * t74 + t122 * t73) * t74) * t73; m(6) * (-t11 * t73 - t12 * t74); 0; m(6) * (-t15 * t73 - t16 * t74); m(6) * (t31 * t74 - t32 * t73); m(6) * t112;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
