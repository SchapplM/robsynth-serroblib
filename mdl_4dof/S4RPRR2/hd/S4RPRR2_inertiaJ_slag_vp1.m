% Calculate joint inertia matrix for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (624->67), mult. (458->95), div. (0->0), fcn. (376->8), ass. (0->45)
t46 = qJ(1) + pkin(7);
t44 = qJ(3) + t46;
t40 = sin(t44);
t72 = t40 ^ 2;
t41 = cos(t44);
t71 = t41 ^ 2;
t47 = sin(qJ(4));
t75 = Icges(5,5) * t47;
t49 = cos(qJ(4));
t74 = Icges(5,6) * t49;
t27 = t74 + t75;
t73 = t40 * t41;
t30 = t47 * rSges(5,1) + t49 * rSges(5,2);
t68 = m(5) * t30;
t48 = sin(qJ(1));
t67 = t48 * pkin(1);
t66 = rSges(5,1) * t49;
t65 = rSges(5,2) * t47;
t64 = t41 * rSges(5,3) + t40 * t65;
t43 = cos(t46);
t50 = cos(qJ(1));
t45 = t50 * pkin(1);
t63 = pkin(2) * t43 + t45;
t60 = Icges(5,2) * t49 ^ 2 + Icges(4,3) + (Icges(5,1) * t47 + 0.2e1 * Icges(5,4) * t49) * t47;
t59 = t27 * t71 + (t75 / 0.2e1 + t74 / 0.2e1 + t27 / 0.2e1) * t72;
t21 = t41 * rSges(4,1) - t40 * rSges(4,2);
t42 = sin(t46);
t58 = -pkin(2) * t42 - t67;
t20 = -t40 * rSges(4,1) - t41 * rSges(4,2);
t52 = Icges(5,5) * t49 - Icges(5,6) * t47;
t51 = t40 * rSges(5,3) + (-t65 + t66) * t41;
t9 = t41 * pkin(3) + t40 * pkin(6) + t51;
t8 = t41 * pkin(6) + (-pkin(3) - t66) * t40 + t64;
t32 = t50 * rSges(2,1) - t48 * rSges(2,2);
t31 = -t48 * rSges(2,1) - t50 * rSges(2,2);
t19 = t43 * rSges(3,1) - t42 * rSges(3,2) + t45;
t18 = -t42 * rSges(3,1) - t43 * rSges(3,2) - t67;
t17 = t21 + t63;
t16 = t20 + t58;
t11 = Icges(5,3) * t40 + t52 * t41;
t10 = -Icges(5,3) * t41 + t52 * t40;
t7 = t9 + t63;
t6 = t58 + t8;
t3 = t40 * (t40 * t66 - t64) + t41 * t51;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t31 ^ 2 + t32 ^ 2) + m(3) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t60; 0; m(3) + m(4) + m(5); m(4) * (t20 * t16 + t21 * t17) + m(5) * (t8 * t6 + t9 * t7) + t60; 0; m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t60; (-t40 * t7 - t41 * t6) * t68 + t59; m(5) * t3; (-t40 * t9 - t41 * t8) * t68 + t59; m(5) * (t3 ^ 2 + (t71 + t72) * t30 ^ 2) + t40 * (-t10 * t73 + t72 * t11) - t41 * (t71 * t10 - t11 * t73);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
