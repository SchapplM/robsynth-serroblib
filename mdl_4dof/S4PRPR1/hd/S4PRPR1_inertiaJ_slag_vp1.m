% Calculate joint inertia matrix for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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

function Mq = S4PRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (186->39), mult. (170->51), div. (0->0), fcn. (144->4), ass. (0->18)
t22 = cos(qJ(4));
t21 = sin(qJ(4));
t19 = pkin(6) + qJ(2);
t17 = sin(t19);
t18 = cos(t19);
t20 = t18 * pkin(2) + t17 * qJ(3);
t7 = -t17 * t21 - t18 * t22;
t8 = t17 * t22 - t18 * t21;
t3 = t8 * rSges(5,1) + t7 * rSges(5,2);
t4 = t7 * rSges(5,1) - t8 * rSges(5,2);
t15 = t18 * qJ(3);
t10 = t18 * rSges(3,1) - t17 * rSges(3,2);
t9 = -t17 * rSges(3,1) - t18 * rSges(3,2);
t6 = t18 * rSges(4,1) + t17 * rSges(4,3) + t20;
t5 = t18 * rSges(4,3) + t15 + (-rSges(4,1) - pkin(2)) * t17;
t2 = t18 * pkin(3) + t20 - t4;
t1 = t15 + (-pkin(2) - pkin(3)) * t17 - t3;
t11 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + Icges(4,2) + Icges(5,3) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; m(4) * (t17 * t5 - t18 * t6) + m(5) * (t17 * t1 - t18 * t2); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t17 ^ 2 + t18 ^ 2); 0; m(5) * (t3 * t1 + t4 * t2) - Icges(5,3); m(5) * (t3 * t17 - t4 * t18); m(5) * (t3 ^ 2 + t4 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
