% Calculate joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:28
% DurationCPUTime: 0.61s
% Computational Cost: add. (620->88), mult. (1296->132), div. (0->0), fcn. (1492->6), ass. (0->45)
t52 = sin(qJ(4));
t53 = cos(qJ(4));
t95 = Icges(5,6) + Icges(6,6);
t96 = Icges(5,5) + Icges(6,5);
t94 = t52 * t96 + t53 * t95;
t98 = t94 / 0.2e1;
t49 = sin(pkin(7));
t50 = cos(pkin(7));
t81 = sin(qJ(3));
t82 = cos(qJ(3));
t32 = -t49 * t81 - t50 * t82;
t31 = t32 ^ 2;
t93 = t52 * t95 - t53 * t96;
t90 = Icges(5,3) + Icges(6,3);
t75 = -rSges(6,3) - qJ(5) - pkin(6);
t33 = -t49 * t82 + t50 * t81;
t87 = t32 * t90 - t33 * t93;
t86 = t32 * t93 + t33 * t90;
t30 = t33 ^ 2;
t43 = -rSges(5,1) * t52 - rSges(5,2) * t53;
t83 = m(5) * t43;
t80 = rSges(5,2) * t52;
t79 = rSges(6,2) * t52;
t78 = t33 * rSges(5,3);
t77 = t33 * t52;
t76 = t33 * t53;
t74 = -rSges(5,1) * t76 - rSges(5,3) * t32;
t73 = t30 + t31;
t68 = t53 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t52;
t48 = pkin(4) * t53 + pkin(3);
t67 = -rSges(6,1) * t76 + t32 * t75 - t33 * t48;
t65 = -rSges(5,1) * t53 + t80;
t54 = rSges(6,1) * t53 + t48 - t79;
t29 = t33 * pkin(6);
t22 = rSges(4,1) * t32 + rSges(4,2) * t33;
t21 = -rSges(4,1) * t33 + rSges(4,2) * t32;
t20 = t68 * t32;
t19 = t68 * t33;
t6 = -t78 - t29 + (pkin(3) - t65) * t32;
t5 = -t32 * pkin(6) + (-pkin(3) + t80) * t33 + t74;
t4 = t54 * t32 + t33 * t75;
t3 = rSges(6,2) * t77 + t67;
t2 = t33 * (rSges(5,2) * t77 + t74) + (t32 * t65 + t78) * t32;
t1 = ((pkin(3) + t79) * t33 + t67) * t33 + (-t29 + (pkin(3) - t54) * t32 + (pkin(6) - t75) * t33) * t32;
t7 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t49 ^ 2 + t50 ^ 2); 0; m(4) * (t21 * t49 - t22 * t50) + m(5) * (t49 * t5 - t50 * t6) + m(6) * (t3 * t49 - t4 * t50); Icges(4,3) + m(4) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (Icges(5,2) + Icges(6,2)) * t53 ^ 2 + (0.2e1 * (Icges(6,4) + Icges(5,4)) * t53 + (Icges(5,1) + Icges(6,1)) * t52) * t52; m(5) * t2 + m(6) * t1; m(6) * (-t19 * t50 + t20 * t49) + (-t32 * t49 + t33 * t50) * t83; m(6) * (t19 * t4 + t20 * t3) - t32 * t5 * t83 + t31 * t98 - t94 * (-t30 / 0.2e1 - t31 / 0.2e1) + (t33 * t98 - t6 * t83) * t33; m(5) * (t43 ^ 2 * t73 + t2 ^ 2) + m(6) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t86 * t33 * t30 + (t87 * t31 + (t32 * t86 + t33 * t87) * t33) * t32; 0; m(6) * (t32 * t50 + t33 * t49); m(6) * (t3 * t33 - t32 * t4); m(6) * (-t19 * t32 + t20 * t33); m(6) * t73;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
