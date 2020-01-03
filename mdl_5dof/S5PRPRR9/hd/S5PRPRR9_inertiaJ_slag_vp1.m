% Calculate joint inertia matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:38
% DurationCPUTime: 0.39s
% Computational Cost: add. (1044->94), mult. (1038->143), div. (0->0), fcn. (1146->6), ass. (0->49)
t47 = sin(qJ(5));
t79 = t47 / 0.2e1;
t48 = cos(qJ(5));
t78 = t48 / 0.2e1;
t77 = Icges(6,5) * t79 + Icges(6,6) * t78;
t46 = pkin(8) + qJ(2);
t44 = sin(t46);
t45 = cos(t46);
t70 = sin(qJ(4));
t71 = cos(qJ(4));
t26 = -t44 * t70 - t45 * t71;
t27 = -t44 * t71 + t45 * t70;
t76 = t27 * t26;
t68 = rSges(6,2) * t47;
t59 = -pkin(4) + t68;
t69 = rSges(6,1) * t48;
t66 = t27 * rSges(6,3) - t26 * t69;
t7 = -t27 * pkin(7) - t59 * t26 - t66;
t67 = t26 * rSges(6,3) + t27 * t69;
t6 = -t26 * pkin(7) + t59 * t27 - t67;
t63 = Icges(6,4) * t48;
t51 = Icges(6,2) * t47 - t63;
t64 = Icges(6,4) * t47;
t52 = -Icges(6,1) * t48 + t64;
t75 = -(t26 * t77 - t48 * (-Icges(6,6) * t26 + t51 * t27) / 0.2e1 - t47 * (-Icges(6,5) * t26 + t52 * t27) / 0.2e1) * t26 - (t27 * t77 + (Icges(6,6) * t27 + t51 * t26) * t78 + (Icges(6,5) * t27 + t52 * t26) * t79) * t27;
t74 = t26 ^ 2;
t73 = t27 ^ 2;
t36 = -t47 * rSges(6,1) - t48 * rSges(6,2);
t72 = m(6) * t36;
t65 = t45 * pkin(2) + t44 * qJ(3);
t60 = t45 * pkin(3) + t65;
t41 = t45 * qJ(3);
t56 = t41 + (-pkin(2) - pkin(3)) * t44;
t18 = -t27 * rSges(5,1) + t26 * rSges(5,2);
t19 = t26 * rSges(5,1) + t27 * rSges(5,2);
t50 = -Icges(6,5) * t48 + Icges(6,6) * t47;
t49 = t48 * (-Icges(6,2) * t48 - t64) + t47 * (-Icges(6,1) * t47 - t63) - Icges(5,3);
t30 = t45 * rSges(3,1) - t44 * rSges(3,2);
t29 = -t44 * rSges(3,1) - t45 * rSges(3,2);
t21 = t45 * rSges(4,1) + t44 * rSges(4,3) + t65;
t20 = t45 * rSges(4,3) + t41 + (-rSges(4,1) - pkin(2)) * t44;
t17 = -t19 + t60;
t16 = -t18 + t56;
t11 = Icges(6,3) * t27 + t50 * t26;
t10 = -Icges(6,3) * t26 + t50 * t27;
t5 = t60 - t7;
t4 = t56 - t6;
t1 = t27 * (t27 * t68 - t67) + t26 * (t26 * t68 + t66);
t2 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(4,2) + Icges(3,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t29 ^ 2 + t30 ^ 2) - t49; 0; m(6) * (t44 * t4 - t45 * t5) + m(5) * (t44 * t16 - t45 * t17) + m(4) * (t44 * t20 - t45 * t21); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t44 ^ 2 + t45 ^ 2); 0; m(6) * (t6 * t4 + t7 * t5) + m(5) * (t18 * t16 + t19 * t17) + t49; m(5) * (t18 * t44 - t19 * t45) + m(6) * (t6 * t44 - t7 * t45); m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) - t49; m(6) * t1; (-t26 * t4 - t27 * t5) * t72 + t75; (-t26 * t44 + t27 * t45) * t72; (-t26 * t6 - t27 * t7) * t72 - t75; m(6) * (t1 ^ 2 + (t73 + t74) * t36 ^ 2) + t27 * (-t10 * t76 + t73 * t11) - t26 * (t74 * t10 - t11 * t76);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
