% Calculate joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:10
% DurationCPUTime: 0.83s
% Computational Cost: add. (1113->123), mult. (2272->188), div. (0->0), fcn. (2588->6), ass. (0->62)
t86 = cos(qJ(4));
t151 = t86 ^ 2;
t126 = sin(qJ(3));
t127 = cos(qJ(3));
t85 = sin(qJ(1));
t87 = cos(qJ(1));
t60 = t87 * t126 - t85 * t127;
t130 = t60 ^ 2;
t59 = -t85 * t126 - t87 * t127;
t57 = t59 ^ 2;
t111 = t57 + t130;
t150 = Icges(6,4) + Icges(5,5);
t149 = -t86 / 0.2e1;
t84 = sin(qJ(4));
t147 = -t150 * t86 + (Icges(5,6) - Icges(6,6)) * t84;
t145 = Icges(6,2) + Icges(5,3);
t143 = rSges(6,1) + pkin(4);
t133 = t145 * t59 - t147 * t60;
t132 = t145 * t60 + t147 * t59;
t88 = Icges(4,3) + (Icges(5,2) + Icges(6,3)) * t151 + ((Icges(5,1) + Icges(6,1)) * t84 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t86) * t84;
t131 = t111 * (-t150 * t84 - 0.2e1 * Icges(6,6) * t149 - (t86 / 0.2e1 - t149) * Icges(5,6));
t72 = -t84 * rSges(5,1) - t86 * rSges(5,2);
t129 = m(5) * t72;
t128 = m(6) * t84;
t125 = t59 * t84;
t124 = t59 * t86;
t123 = t60 * t84;
t122 = t60 * t86;
t121 = t60 * pkin(3) + t59 * pkin(7);
t120 = t59 * pkin(3) - t60 * pkin(7);
t118 = (-rSges(6,3) - qJ(5)) * t86 + t143 * t84;
t117 = t87 * pkin(1) + t85 * qJ(2);
t112 = qJ(5) * t84;
t110 = t87 * pkin(2) + t117;
t80 = t87 * qJ(2);
t108 = t80 + (-pkin(1) - pkin(2)) * t85;
t107 = -t59 * rSges(6,2) - rSges(6,3) * t123 - t60 * t112 - t143 * t122;
t106 = t60 * rSges(6,2) - rSges(6,3) * t125 - t59 * t112 - t143 * t124;
t35 = -t60 * rSges(4,1) + t59 * rSges(4,2);
t36 = t59 * rSges(4,1) + t60 * rSges(4,2);
t101 = -t59 * t85 + t60 * t87;
t92 = -rSges(5,1) * t122 + rSges(5,2) * t123 - t59 * rSges(5,3);
t91 = -rSges(5,1) * t124 + rSges(5,2) * t125 + t60 * rSges(5,3);
t5 = t107 - t121;
t6 = -t106 + t120;
t13 = t92 - t121;
t14 = -t91 + t120;
t74 = t87 * rSges(2,1) - t85 * rSges(2,2);
t73 = -t85 * rSges(2,1) - t87 * rSges(2,2);
t38 = t87 * rSges(3,1) + t85 * rSges(3,3) + t117;
t37 = t87 * rSges(3,3) + t80 + (-rSges(3,1) - pkin(1)) * t85;
t34 = -t36 + t110;
t33 = -t35 + t108;
t32 = t118 * t59;
t31 = t118 * t60;
t12 = -t14 + t110;
t11 = -t13 + t108;
t4 = -t6 + t110;
t3 = -t5 + t108;
t2 = t59 * t91 + t60 * t92;
t1 = t106 * t59 + t107 * t60;
t7 = [Icges(3,2) + Icges(2,3) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(2) * (t73 ^ 2 + t74 ^ 2) + t88; m(6) * (t85 * t3 - t87 * t4) + m(5) * (t85 * t11 - t87 * t12) + m(4) * (t85 * t33 - t87 * t34) + m(3) * (t85 * t37 - t87 * t38); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t85 ^ 2 + t87 ^ 2); m(6) * (t5 * t3 + t6 * t4) + m(5) * (t13 * t11 + t14 * t12) + m(4) * (t35 * t33 + t36 * t34) - t88; m(4) * (t35 * t85 - t36 * t87) + m(5) * (t13 * t85 - t14 * t87) + m(6) * (t5 * t85 - t6 * t87); m(6) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + t88; m(6) * (t32 * t3 + t31 * t4) + (-t11 * t59 - t12 * t60) * t129 + t131; m(6) * (-t31 * t87 + t32 * t85) + t101 * t129; m(6) * (t31 * t6 + t32 * t5) + (-t13 * t59 - t14 * t60) * t129 - t131; m(5) * (t111 * t72 ^ 2 + t2 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t132 * t60 * t130 + (t133 * t57 + (t132 * t59 + t133 * t60) * t60) * t59; (-t3 * t59 - t4 * t60) * t128; t101 * t128; (-t5 * t59 - t6 * t60) * t128; m(6) * (t86 * t1 + (-t31 * t60 - t32 * t59) * t84); m(6) * (t111 * t84 ^ 2 + t151);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
