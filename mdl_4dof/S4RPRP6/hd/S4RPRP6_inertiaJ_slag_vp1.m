% Calculate joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.51s
% Computational Cost: add. (301->79), mult. (652->114), div. (0->0), fcn. (544->4), ass. (0->40)
t50 = sin(qJ(1));
t46 = t50 ^ 2;
t52 = cos(qJ(1));
t47 = t52 ^ 2;
t35 = t46 + t47;
t91 = Icges(4,5) + Icges(5,5);
t90 = Icges(4,6) + Icges(5,6);
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t82 = rSges(5,1) + pkin(3);
t89 = (rSges(5,2) * t51 + t82 * t49) * t52;
t88 = t91 * t49 + t90 * t51;
t85 = Icges(4,3) + Icges(5,3);
t81 = t88 * t50 + t85 * t52;
t80 = (rSges(4,1) * t49 + rSges(4,2) * t51) * t52;
t79 = t85 * t50 - t88 * t52;
t75 = t49 * t50;
t74 = t50 * t51;
t73 = t52 * pkin(1) + t50 * qJ(2);
t68 = rSges(4,1) * t75 + rSges(4,2) * t74 + t52 * rSges(4,3);
t67 = -t49 * rSges(5,2) + t82 * t51;
t66 = -rSges(5,2) * t74 - t52 * rSges(5,3) - t82 * t75;
t42 = t52 * qJ(2);
t5 = t42 + t80 + (-rSges(4,3) - pkin(1) - pkin(5)) * t50;
t6 = t52 * pkin(5) + t68 + t73;
t53 = m(4) * (t50 * t5 - t52 * t6);
t48 = -qJ(4) - pkin(5);
t34 = t52 * rSges(2,1) - t50 * rSges(2,2);
t33 = t51 * rSges(4,1) - t49 * rSges(4,2);
t31 = -t50 * rSges(2,1) - t52 * rSges(2,2);
t23 = m(5) * t35;
t22 = -t52 * rSges(3,2) + t50 * rSges(3,3) + t73;
t21 = t52 * rSges(3,3) + t42 + (rSges(3,2) - pkin(1)) * t50;
t8 = t67 * t52;
t7 = t67 * t50;
t4 = -t52 * t48 - t66 + t73;
t3 = t42 + t89 + (-rSges(5,3) - pkin(1) + t48) * t50;
t2 = -t50 * t68 + (t50 * rSges(4,3) - t80) * t52;
t1 = t66 * t50 + (t50 * rSges(5,3) - t89) * t52;
t9 = [Icges(3,1) + Icges(2,3) + m(2) * (t31 ^ 2 + t34 ^ 2) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,1) + Icges(5,1)) * t51 ^ 2 + (0.2e1 * (-Icges(5,4) - Icges(4,4)) * t51 + (Icges(4,2) + Icges(5,2)) * t49) * t49; m(3) * (t50 * t21 - t52 * t22) + t53 + m(5) * (t50 * t3 - t52 * t4); t23 + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t35; m(5) * (t7 * t3 - t8 * t4) + t33 * t53 + t35 * (-t90 * t49 + t91 * t51); m(5) * (t7 * t50 + t8 * t52) + m(4) * t35 * t33; m(4) * (t35 * t33 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t79 * t50 * t46 + (t81 * t47 + (t81 * t50 + t79 * t52) * t50) * t52; m(5) * (t52 * t3 + t50 * t4); 0; m(5) * (-t50 * t8 + t52 * t7); t23;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
