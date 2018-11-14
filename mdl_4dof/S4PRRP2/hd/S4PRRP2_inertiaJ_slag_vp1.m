% Calculate joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (117->34), mult. (132->40), div. (0->0), fcn. (66->4), ass. (0->20)
t22 = rSges(5,1) + pkin(3);
t18 = sin(qJ(2));
t21 = pkin(2) * t18;
t20 = Icges(4,3) + Icges(5,3);
t17 = qJ(2) + qJ(3);
t14 = sin(t17);
t15 = cos(t17);
t8 = t15 * rSges(4,1) - rSges(4,2) * t14;
t4 = -rSges(5,2) * t14 + t22 * t15;
t7 = -rSges(4,1) * t14 - rSges(4,2) * t15;
t3 = -rSges(5,2) * t15 - t22 * t14;
t19 = cos(qJ(2));
t16 = t19 * pkin(2);
t10 = rSges(3,1) * t19 - t18 * rSges(3,2);
t9 = -t18 * rSges(3,1) - rSges(3,2) * t19;
t6 = t16 + t8;
t5 = t7 - t21;
t2 = t16 + t4;
t1 = t3 - t21;
t11 = [m(2) + m(3) + m(4) + m(5); m(3) * t10 + m(4) * t6 + m(5) * t2; Icges(3,3) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t20; m(4) * t8 + m(5) * t4; m(4) * (t5 * t7 + t6 * t8) + m(5) * (t1 * t3 + t2 * t4) + t20; m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t20; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
