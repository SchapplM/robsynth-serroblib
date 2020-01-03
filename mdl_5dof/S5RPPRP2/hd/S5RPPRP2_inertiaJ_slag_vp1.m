% Calculate joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:07
% DurationCPUTime: 0.66s
% Computational Cost: add. (994->103), mult. (791->144), div. (0->0), fcn. (669->8), ass. (0->53)
t55 = pkin(8) + qJ(4);
t52 = cos(t55);
t106 = t52 ^ 2;
t56 = qJ(1) + pkin(7);
t51 = sin(t56);
t48 = t51 ^ 2;
t53 = cos(t56);
t49 = t53 ^ 2;
t84 = t48 + t49;
t104 = Icges(6,4) + Icges(5,5);
t103 = Icges(5,6) - Icges(6,6);
t50 = sin(t55);
t101 = -t103 * t50 + t104 * t52;
t99 = Icges(6,2) + Icges(5,3);
t93 = rSges(6,1) + pkin(4);
t95 = t93 * t52;
t94 = rSges(6,3) + qJ(5);
t92 = -t101 * t51 + t99 * t53;
t91 = t101 * t53 + t99 * t51;
t60 = sin(qJ(1));
t88 = t60 * pkin(1);
t87 = t50 * t53;
t86 = t52 * t53;
t85 = -t93 * t50 + t94 * t52;
t79 = qJ(5) * t50;
t78 = rSges(4,3) + qJ(3);
t58 = cos(pkin(8));
t46 = t58 * pkin(3) + pkin(2);
t61 = cos(qJ(1));
t54 = t61 * pkin(1);
t59 = -pkin(6) - qJ(3);
t76 = t53 * t46 - t51 * t59 + t54;
t75 = t51 * rSges(6,2) + rSges(6,3) * t87 + t53 * t79 + t93 * t86;
t74 = rSges(5,1) * t52 - rSges(5,2) * t50;
t63 = rSges(5,1) * t86 - rSges(5,2) * t87 + t51 * rSges(5,3);
t57 = sin(pkin(8));
t62 = rSges(4,1) * t58 - rSges(4,2) * t57 + pkin(2);
t43 = t61 * rSges(2,1) - t60 * rSges(2,2);
t42 = -t60 * rSges(2,1) - t61 * rSges(2,2);
t34 = t50 * rSges(5,1) + t52 * rSges(5,2);
t24 = t53 * rSges(3,1) - t51 * rSges(3,2) + t54;
t23 = -t51 * rSges(3,1) - t53 * rSges(3,2) - t88;
t10 = t85 * t53;
t9 = t85 * t51;
t8 = t78 * t51 + t62 * t53 + t54;
t7 = -t62 * t51 + t78 * t53 - t88;
t6 = t63 + t76;
t5 = -t88 + (rSges(5,3) - t59) * t53 + (-t46 - t74) * t51;
t4 = t53 * t63 + (-t53 * rSges(5,3) + t74 * t51) * t51;
t3 = t75 + t76;
t2 = -t88 + (rSges(6,2) - t59) * t53 + (-t94 * t50 - t46 - t95) * t51;
t1 = t75 * t53 + (-t53 * rSges(6,2) + (rSges(6,3) * t50 + t79 + t95) * t51) * t51;
t11 = [Icges(4,2) * t58 ^ 2 + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t57 + 0.2e1 * Icges(4,4) * t58) * t57 + m(2) * (t42 ^ 2 + t43 ^ 2) + m(3) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,2) + Icges(6,3)) * t106 + (0.2e1 * (-Icges(6,5) + Icges(5,4)) * t52 + (Icges(5,1) + Icges(6,1)) * t50) * t50; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t51 * t7 - t53 * t8) + m(5) * (t51 * t5 - t53 * t6) + m(6) * (t51 * t2 - t53 * t3); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t84; m(6) * (t10 * t2 + t9 * t3) + m(5) * (-t5 * t53 - t51 * t6) * t34 + t84 * (t103 * t52 + t104 * t50); m(5) * t4 + m(6) * t1; m(6) * (t10 * t51 - t9 * t53); m(5) * (t84 * t34 ^ 2 + t4 ^ 2) + m(6) * (t1 ^ 2 + t10 ^ 2 + t9 ^ 2) + t91 * t51 * t48 + (t92 * t49 + (t92 * t51 + t91 * t53) * t51) * t53; m(6) * (t2 * t53 + t3 * t51) * t50; -m(6) * t52; 0; m(6) * (-t52 * t1 + (t10 * t53 + t51 * t9) * t50); m(6) * (t84 * t50 ^ 2 + t106);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
