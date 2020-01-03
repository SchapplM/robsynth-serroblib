% Calculate joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:33
% DurationCPUTime: 0.51s
% Computational Cost: add. (559->82), mult. (584->106), div. (0->0), fcn. (500->6), ass. (0->45)
t46 = qJ(1) + pkin(6);
t43 = sin(t46);
t41 = t43 ^ 2;
t44 = cos(t46);
t42 = t44 ^ 2;
t71 = t41 + t42;
t91 = Icges(4,5) + Icges(5,5);
t90 = Icges(4,6) + Icges(5,6);
t48 = sin(qJ(3));
t50 = cos(qJ(3));
t88 = -t90 * t48 + t91 * t50;
t85 = Icges(4,3) + Icges(5,3);
t82 = -t88 * t43 + t85 * t44;
t81 = t85 * t43 + t88 * t44;
t49 = sin(qJ(1));
t78 = t49 * pkin(1);
t77 = rSges(4,2) * t48;
t76 = rSges(5,2) * t48;
t75 = t44 * rSges(4,3);
t74 = t44 * t48;
t73 = t44 * t50;
t72 = rSges(4,1) * t73 + t43 * rSges(4,3);
t40 = t50 * pkin(3) + pkin(2);
t66 = rSges(5,1) * t73 + t43 * rSges(5,3) + t44 * t40;
t65 = -t50 * rSges(5,2) + (-rSges(5,1) - pkin(3)) * t48;
t63 = rSges(4,1) * t50 - t77;
t52 = -rSges(5,1) * t50 - t40 + t76;
t51 = cos(qJ(1));
t47 = -qJ(4) - pkin(5);
t45 = t51 * pkin(1);
t39 = t44 * pkin(5);
t34 = t51 * rSges(2,1) - t49 * rSges(2,2);
t33 = -t49 * rSges(2,1) - t51 * rSges(2,2);
t32 = t48 * rSges(4,1) + t50 * rSges(4,2);
t22 = t44 * rSges(3,1) - t43 * rSges(3,2) + t45;
t21 = -t43 * rSges(3,1) - t44 * rSges(3,2) - t78;
t20 = t65 * t44;
t19 = t65 * t43;
t6 = t43 * pkin(5) + t45 + (pkin(2) - t77) * t44 + t72;
t5 = t75 - t78 + t39 + (-pkin(2) - t63) * t43;
t4 = -rSges(5,2) * t74 - t43 * t47 + t45 + t66;
t3 = -t78 + (rSges(5,3) - t47) * t44 + t52 * t43;
t2 = t44 * (-rSges(4,2) * t74 + t72) + (t63 * t43 - t75) * t43;
t1 = ((-pkin(2) - t76) * t44 + t66) * t44 + (t39 + (-pkin(2) - t52) * t43 + (-pkin(5) - rSges(5,3)) * t44) * t43;
t7 = [Icges(2,3) + Icges(3,3) + m(2) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,2) + Icges(5,2)) * t50 ^ 2 + (0.2e1 * (Icges(5,4) + Icges(4,4)) * t50 + (Icges(4,1) + Icges(5,1)) * t48) * t48; 0; m(3) + m(4) + m(5); m(5) * (t19 * t4 + t20 * t3) + m(4) * (-t43 * t6 - t44 * t5) * t32 + t71 * (t91 * t48 + t90 * t50); m(4) * t2 + m(5) * t1; m(4) * (t71 * t32 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t81 * t43 * t41 + (t82 * t42 + (t82 * t43 + t81 * t44) * t43) * t44; m(5) * (t43 * t3 - t44 * t4); 0; m(5) * (-t44 * t19 + t43 * t20); m(5) * t71;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
