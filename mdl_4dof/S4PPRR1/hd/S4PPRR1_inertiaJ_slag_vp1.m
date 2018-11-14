% Calculate joint inertia matrix for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:08
% EndTime: 2018-11-14 13:40:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (141->32), mult. (186->49), div. (0->0), fcn. (170->6), ass. (0->21)
t18 = sin(pkin(6));
t19 = cos(pkin(6));
t17 = qJ(3) + qJ(4);
t25 = sin(t17);
t26 = cos(t17);
t10 = t18 * t26 - t19 * t25;
t9 = -t18 * t25 - t19 * t26;
t3 = t10 * rSges(5,1) + t9 * rSges(5,2);
t21 = cos(qJ(3));
t24 = t21 * pkin(3);
t20 = sin(qJ(3));
t23 = t18 * t20;
t22 = t19 * t20;
t4 = t9 * rSges(5,1) - t10 * rSges(5,2);
t12 = t18 * t21 - t22;
t11 = -t19 * t21 - t23;
t6 = t11 * rSges(4,1) - t12 * rSges(4,2);
t5 = t12 * rSges(4,1) + t11 * rSges(4,2);
t2 = -pkin(3) * t23 - t24 * t19 + t4;
t1 = -pkin(3) * t22 + t24 * t18 + t3;
t7 = [m(2) + m(3) + m(4) + m(5); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t18 ^ 2 + t19 ^ 2); 0; m(4) * (t5 * t18 - t6 * t19) + m(5) * (t1 * t18 - t2 * t19); Icges(4,3) + Icges(5,3) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; m(5) * (t3 * t18 - t4 * t19); m(5) * (t3 * t1 + t4 * t2) + Icges(5,3); m(5) * (t3 ^ 2 + t4 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
