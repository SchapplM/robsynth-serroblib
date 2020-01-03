% Calculate joint inertia matrix for
% S5RPPRP3
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (702->93), mult. (726->125), div. (0->0), fcn. (602->6), ass. (0->47)
t53 = qJ(1) + pkin(7);
t50 = sin(t53);
t48 = t50 ^ 2;
t51 = cos(t53);
t49 = t51 ^ 2;
t26 = t48 + t49;
t99 = Icges(5,5) + Icges(6,5);
t98 = Icges(5,6) + Icges(6,6);
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t90 = rSges(6,1) + pkin(4);
t97 = (rSges(6,2) * t57 + t90 * t55) * t51;
t96 = t99 * t55 + t98 * t57;
t93 = Icges(5,3) + Icges(6,3);
t89 = t96 * t50 + t93 * t51;
t88 = (rSges(5,1) * t55 + rSges(5,2) * t57) * t51;
t87 = t93 * t50 - t96 * t51;
t56 = sin(qJ(1));
t84 = t56 * pkin(1);
t82 = t50 * t55;
t81 = t50 * t57;
t76 = rSges(5,1) * t82 + rSges(5,2) * t81 + t51 * rSges(5,3);
t58 = cos(qJ(1));
t52 = t58 * pkin(1);
t75 = t51 * pkin(2) + t50 * qJ(3) + t52;
t74 = t51 * qJ(3) - t84;
t73 = -t55 * rSges(6,2) + t90 * t57;
t72 = -rSges(6,2) * t81 - t51 * rSges(6,3) - t90 * t82;
t5 = t88 + (-rSges(5,3) - pkin(2) - pkin(6)) * t50 + t74;
t6 = t51 * pkin(6) + t75 + t76;
t59 = m(5) * (t50 * t5 - t51 * t6);
t54 = -qJ(5) - pkin(6);
t37 = t58 * rSges(2,1) - t56 * rSges(2,2);
t36 = t57 * rSges(5,1) - t55 * rSges(5,2);
t34 = -t56 * rSges(2,1) - t58 * rSges(2,2);
t25 = m(6) * t26;
t24 = t51 * rSges(3,1) - t50 * rSges(3,2) + t52;
t23 = -t50 * rSges(3,1) - t51 * rSges(3,2) - t84;
t22 = t73 * t51;
t21 = t73 * t50;
t8 = -t51 * rSges(4,2) + t50 * rSges(4,3) + t75;
t7 = t51 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t50 + t74;
t4 = -t51 * t54 - t72 + t75;
t3 = t97 + (-rSges(6,3) - pkin(2) + t54) * t50 + t74;
t2 = -t50 * t76 + (t50 * rSges(5,3) - t88) * t51;
t1 = t72 * t50 + (t50 * rSges(6,3) - t97) * t51;
t9 = [Icges(4,1) + Icges(2,3) + Icges(3,3) + m(2) * (t34 ^ 2 + t37 ^ 2) + m(3) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (Icges(5,1) + Icges(6,1)) * t57 ^ 2 + (0.2e1 * (-Icges(6,4) - Icges(5,4)) * t57 + (Icges(5,2) + Icges(6,2)) * t55) * t55; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t50 * t7 - t51 * t8) + t59 + m(6) * (t50 * t3 - t51 * t4); 0; t25 + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t26; m(6) * (t21 * t3 - t22 * t4) + t36 * t59 + t26 * (-t98 * t55 + t99 * t57); m(5) * t2 + m(6) * t1; m(6) * (t21 * t50 + t22 * t51) + m(5) * t26 * t36; m(5) * (t26 * t36 ^ 2 + t2 ^ 2) + m(6) * (t1 ^ 2 + t21 ^ 2 + t22 ^ 2) + t87 * t50 * t48 + (t89 * t49 + (t89 * t50 + t87 * t51) * t50) * t51; m(6) * (t51 * t3 + t50 * t4); 0; 0; m(6) * (t51 * t21 - t50 * t22); t25;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
