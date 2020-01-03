% Calculate joint inertia matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:34
% DurationCPUTime: 0.41s
% Computational Cost: add. (495->89), mult. (1006->142), div. (0->0), fcn. (1110->6), ass. (0->48)
t44 = sin(qJ(4));
t78 = t44 / 0.2e1;
t46 = cos(qJ(4));
t77 = t46 / 0.2e1;
t76 = Icges(5,5) * t78 + Icges(5,6) * t77;
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t69 = sin(qJ(3));
t70 = cos(qJ(3));
t26 = -t45 * t69 - t47 * t70;
t27 = -t45 * t70 + t47 * t69;
t75 = t27 * t26;
t67 = rSges(5,2) * t44;
t58 = -pkin(3) + t67;
t68 = rSges(5,1) * t46;
t65 = t27 * rSges(5,3) - t26 * t68;
t7 = -t27 * pkin(6) - t58 * t26 - t65;
t66 = t26 * rSges(5,3) + t27 * t68;
t6 = -t26 * pkin(6) + t58 * t27 - t66;
t62 = Icges(5,4) * t46;
t50 = Icges(5,2) * t44 - t62;
t63 = Icges(5,4) * t44;
t51 = -Icges(5,1) * t46 + t63;
t74 = -(t26 * t76 - t46 * (-Icges(5,6) * t26 + t50 * t27) / 0.2e1 - t44 * (-Icges(5,5) * t26 + t51 * t27) / 0.2e1) * t26 - (t27 * t76 + (Icges(5,6) * t27 + t50 * t26) * t77 + (Icges(5,5) * t27 + t51 * t26) * t78) * t27;
t73 = t26 ^ 2;
t72 = t27 ^ 2;
t33 = -t44 * rSges(5,1) - t46 * rSges(5,2);
t71 = m(5) * t33;
t64 = t47 * pkin(1) + t45 * qJ(2);
t59 = t47 * pkin(2) + t64;
t41 = t47 * qJ(2);
t55 = t41 + (-pkin(1) - pkin(2)) * t45;
t18 = -t27 * rSges(4,1) + t26 * rSges(4,2);
t19 = t26 * rSges(4,1) + t27 * rSges(4,2);
t49 = -Icges(5,5) * t46 + Icges(5,6) * t44;
t48 = t46 * (-Icges(5,2) * t46 - t63) + t44 * (-Icges(5,1) * t44 - t62) - Icges(4,3);
t35 = t47 * rSges(2,1) - t45 * rSges(2,2);
t34 = -t45 * rSges(2,1) - t47 * rSges(2,2);
t21 = t47 * rSges(3,1) + t45 * rSges(3,3) + t64;
t20 = t47 * rSges(3,3) + t41 + (-rSges(3,1) - pkin(1)) * t45;
t17 = -t19 + t59;
t16 = -t18 + t55;
t11 = Icges(5,3) * t27 + t49 * t26;
t10 = -Icges(5,3) * t26 + t49 * t27;
t5 = t59 - t7;
t4 = t55 - t6;
t1 = t27 * (t27 * t67 - t66) + t26 * (t26 * t67 + t65);
t2 = [Icges(3,2) + Icges(2,3) + m(2) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) - t48; m(3) * (t45 * t20 - t47 * t21) + m(4) * (t45 * t16 - t47 * t17) + m(5) * (t45 * t4 - t47 * t5); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t45 ^ 2 + t47 ^ 2); m(4) * (t18 * t16 + t19 * t17) + m(5) * (t6 * t4 + t7 * t5) + t48; m(4) * (t18 * t45 - t19 * t47) + m(5) * (t6 * t45 - t7 * t47); m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) - t48; (-t26 * t4 - t27 * t5) * t71 + t74; (-t26 * t45 + t27 * t47) * t71; (-t26 * t6 - t27 * t7) * t71 - t74; m(5) * (t1 ^ 2 + (t72 + t73) * t33 ^ 2) + t27 * (-t10 * t75 + t72 * t11) - t26 * (t73 * t10 - t11 * t75);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
