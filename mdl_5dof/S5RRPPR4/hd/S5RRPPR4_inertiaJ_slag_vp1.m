% Calculate joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:42
% DurationCPUTime: 0.44s
% Computational Cost: add. (1244->106), mult. (1188->149), div. (0->0), fcn. (1262->8), ass. (0->54)
t60 = qJ(1) + qJ(2);
t57 = sin(t60);
t58 = cos(t60);
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t33 = -t57 * t77 - t58 * t78;
t89 = t33 ^ 2;
t34 = -t57 * t78 + t58 * t77;
t88 = t34 ^ 2;
t61 = sin(qJ(5));
t92 = Icges(6,5) * t61;
t63 = cos(qJ(5));
t91 = Icges(6,6) * t63;
t40 = -t91 - t92;
t90 = t34 * t33;
t43 = -t61 * rSges(6,1) - t63 * rSges(6,2);
t85 = m(6) * t43;
t62 = sin(qJ(1));
t84 = t62 * pkin(1);
t83 = rSges(6,1) * t63;
t82 = rSges(6,2) * t61;
t81 = t58 * pkin(2) + t57 * qJ(3);
t76 = t58 * pkin(3) + t81;
t75 = t40 * t88 + (-t92 / 0.2e1 - t91 / 0.2e1 + t40 / 0.2e1) * t89;
t37 = t58 * rSges(3,1) - t57 * rSges(3,2);
t74 = t33 * rSges(6,3) - t34 * t82;
t23 = t58 * rSges(4,1) + t57 * rSges(4,3) + t81;
t50 = t58 * qJ(3);
t73 = t50 + (-pkin(2) - pkin(3)) * t57;
t19 = -t33 * rSges(5,1) - t34 * rSges(5,2) + t76;
t36 = -t57 * rSges(3,1) - t58 * rSges(3,2);
t67 = -Icges(6,5) * t63 + Icges(6,6) * t61;
t66 = t34 * rSges(6,3) + (t82 - t83) * t33;
t22 = t50 + t58 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t57;
t65 = t63 ^ 2 * Icges(6,2) + Icges(4,2) + Icges(3,3) + Icges(5,3) + (Icges(6,1) * t61 + 0.2e1 * Icges(6,4) * t63) * t61;
t18 = t34 * rSges(5,1) - t33 * rSges(5,2) + t73;
t7 = -t33 * pkin(4) + t34 * pkin(7) + t66 + t76;
t6 = t33 * pkin(7) + (pkin(4) + t83) * t34 + t73 + t74;
t64 = cos(qJ(1));
t59 = t64 * pkin(1);
t45 = t64 * rSges(2,1) - t62 * rSges(2,2);
t44 = -t62 * rSges(2,1) - t64 * rSges(2,2);
t32 = t37 + t59;
t31 = t36 - t84;
t21 = t59 + t23;
t20 = t22 - t84;
t17 = t59 + t19;
t16 = t18 - t84;
t11 = Icges(6,3) * t34 + t67 * t33;
t10 = -Icges(6,3) * t33 + t67 * t34;
t5 = t59 + t7;
t4 = t6 - t84;
t1 = t34 * (-t34 * t83 - t74) + t33 * t66;
t2 = [Icges(2,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(3) * (t31 ^ 2 + t32 ^ 2) + m(2) * (t44 ^ 2 + t45 ^ 2) + t65; m(6) * (t6 * t4 + t7 * t5) + m(4) * (t22 * t20 + t23 * t21) + m(5) * (t18 * t16 + t19 * t17) + m(3) * (t36 * t31 + t37 * t32) + t65; m(6) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + t65; m(6) * (t57 * t4 - t58 * t5) + m(4) * (t57 * t20 - t58 * t21) + m(5) * (t57 * t16 - t58 * t17); m(6) * (t57 * t6 - t58 * t7) + m(4) * (t57 * t22 - t58 * t23) + m(5) * (t57 * t18 - t58 * t19); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t57 ^ 2 + t58 ^ 2); 0; 0; 0; m(5) + m(6); (-t33 * t4 - t34 * t5) * t85 + t75; (-t33 * t6 - t34 * t7) * t85 + t75; (-t33 * t57 + t34 * t58) * t85; -m(6) * t1; m(6) * (t1 ^ 2 + (t88 + t89) * t43 ^ 2) + t34 * (-t10 * t90 + t88 * t11) - t33 * (t89 * t10 - t11 * t90);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
