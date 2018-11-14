% Calculate joint inertia matrix for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:27
% EndTime: 2018-11-14 13:43:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (234->36), mult. (138->44), div. (0->0), fcn. (84->4), ass. (0->22)
t28 = rSges(5,1) + pkin(3);
t27 = qJ(4) + rSges(5,3);
t24 = pkin(6) + qJ(2);
t21 = sin(t24);
t26 = pkin(2) * t21;
t25 = Icges(4,3) + Icges(5,2);
t23 = qJ(3) + t24;
t19 = sin(t23);
t20 = cos(t23);
t8 = t20 * rSges(4,1) - t19 * rSges(4,2);
t4 = t27 * t19 + t28 * t20;
t7 = -t19 * rSges(4,1) - t20 * rSges(4,2);
t3 = -t28 * t19 + t27 * t20;
t22 = cos(t24);
t18 = pkin(2) * t22;
t10 = t22 * rSges(3,1) - t21 * rSges(3,2);
t9 = -t21 * rSges(3,1) - t22 * rSges(3,2);
t6 = t18 + t8;
t5 = t7 - t26;
t2 = t18 + t4;
t1 = t3 - t26;
t11 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t25; 0; m(4) * (t7 * t5 + t8 * t6) + m(5) * (t3 * t1 + t4 * t2) + t25; m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t25; 0; m(5) * (t19 * t1 - t20 * t2); m(5) * (t19 * t3 - t20 * t4); m(5) * (t19 ^ 2 + t20 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
