% Calculate joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:16
% EndTime: 2019-12-31 17:32:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (485->70), mult. (806->116), div. (0->0), fcn. (900->8), ass. (0->41)
t35 = pkin(8) + qJ(5);
t33 = sin(t35);
t63 = Icges(6,5) * t33;
t34 = cos(t35);
t62 = Icges(6,6) * t34;
t61 = t63 / 0.2e1 + t62 / 0.2e1;
t48 = m(5) / 0.2e1 + m(6) / 0.2e1;
t60 = 0.2e1 * t48;
t37 = sin(pkin(7));
t39 = cos(pkin(7));
t56 = sin(qJ(3));
t57 = cos(qJ(3));
t25 = -t37 * t56 - t39 * t57;
t26 = -t37 * t57 + t39 * t56;
t59 = t26 * t25;
t24 = t25 ^ 2;
t23 = t26 ^ 2;
t22 = -t33 * rSges(6,1) - t34 * rSges(6,2);
t58 = m(6) * t22;
t55 = rSges(6,1) * t34;
t54 = rSges(6,2) * t33;
t53 = -t25 * rSges(6,3) - t26 * t55;
t52 = t23 + t24;
t49 = -rSges(5,3) - qJ(4);
t45 = t54 - t55;
t42 = -Icges(6,5) * t34 + Icges(6,6) * t33;
t36 = sin(pkin(8));
t38 = cos(pkin(8));
t41 = rSges(5,1) * t38 - rSges(5,2) * t36 + pkin(3);
t40 = -pkin(6) - qJ(4);
t32 = t38 * pkin(4) + pkin(3);
t14 = t25 * rSges(4,1) + t26 * rSges(4,2);
t13 = -t26 * rSges(4,1) + t25 * rSges(4,2);
t7 = Icges(6,3) * t26 + t42 * t25;
t6 = -Icges(6,3) * t25 + t42 * t26;
t5 = t41 * t25 + t49 * t26;
t4 = t49 * t25 - t41 * t26;
t3 = (-rSges(6,3) + t40) * t26 + (t32 - t45) * t25;
t2 = t25 * t40 + (-t32 + t54) * t26 + t53;
t1 = t26 * (t26 * t54 + t53) + (t26 * rSges(6,3) + t45 * t25) * t25;
t8 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t48) * (t37 ^ 2 + t39 ^ 2); 0; m(4) * (t13 * t37 - t14 * t39) + m(5) * (t4 * t37 - t5 * t39) + m(6) * (t2 * t37 - t3 * t39); t38 ^ 2 * Icges(5,2) + t34 ^ 2 * Icges(6,2) + Icges(4,3) + (Icges(5,1) * t36 + 0.2e1 * Icges(5,4) * t38) * t36 + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t13 ^ 2 + t14 ^ 2) + (Icges(6,1) * t33 + 0.2e1 * Icges(6,4) * t34) * t33; 0; (t25 * t39 + t26 * t37) * t60; m(6) * (t26 * t2 - t25 * t3) + m(5) * (-t25 * t5 + t26 * t4); t52 * t60; m(6) * t1; (-t25 * t37 + t26 * t39) * t58; (-t24 / 0.2e1 - t23 / 0.2e1) * (-t62 - t63) + (t61 * t26 - t3 * t58) * t26 + (-t2 * t58 + t61 * t25) * t25; 0; m(6) * (t52 * t22 ^ 2 + t1 ^ 2) + t26 * (t23 * t7 - t6 * t59) - t25 * (t24 * t6 - t7 * t59);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
