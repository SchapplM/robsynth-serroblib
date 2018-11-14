% Calculate joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (179->55), mult. (242->71), div. (0->0), fcn. (210->6), ass. (0->27)
t22 = pkin(6) + qJ(4);
t31 = cos(t22);
t30 = sin(t22);
t23 = sin(pkin(6));
t25 = sin(qJ(1));
t29 = t25 * t23;
t26 = cos(qJ(1));
t28 = t26 * t23;
t27 = t26 * pkin(1) + t25 * qJ(2);
t10 = t25 * t31 - t26 * t30;
t9 = -t25 * t30 - t26 * t31;
t6 = t9 * rSges(5,1) - t10 * rSges(5,2);
t5 = t10 * rSges(5,1) + t9 * rSges(5,2);
t24 = cos(pkin(6));
t20 = t26 * qJ(2);
t18 = t24 * pkin(3) + pkin(2);
t14 = t26 * rSges(2,1) - t25 * rSges(2,2);
t13 = -t25 * rSges(2,1) - t26 * rSges(2,2);
t12 = t25 * t24 - t28;
t11 = -t26 * t24 - t29;
t8 = t26 * rSges(3,1) + t25 * rSges(3,3) + t27;
t7 = t26 * rSges(3,3) + t20 + (-rSges(3,1) - pkin(1)) * t25;
t4 = -t11 * rSges(4,1) + t12 * rSges(4,2) + t26 * pkin(2) + t27;
t3 = -t12 * rSges(4,1) - t11 * rSges(4,2) + t20 + (-pkin(1) - pkin(2)) * t25;
t2 = pkin(3) * t29 + t26 * t18 + t27 - t6;
t1 = pkin(3) * t28 + t20 + (-pkin(1) - t18) * t25 - t5;
t15 = [Icges(2,3) + Icges(3,2) + Icges(4,3) + Icges(5,3) + m(2) * (t13 ^ 2 + t14 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); m(3) * (t25 * t7 - t26 * t8) + m(4) * (t25 * t3 - t26 * t4) + m(5) * (t25 * t1 - t26 * t2); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t25 ^ 2 + t26 ^ 2); 0; 0; m(4) + m(5); -Icges(5,3) + m(5) * (t5 * t1 + t6 * t2); m(5) * (t5 * t25 - t6 * t26); 0; m(5) * (t5 ^ 2 + t6 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;
