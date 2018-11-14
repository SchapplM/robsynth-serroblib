% Calculate joint inertia matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:29
% EndTime: 2018-11-14 13:51:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (295->48), mult. (196->59), div. (0->0), fcn. (124->6), ass. (0->30)
t37 = rSges(5,1) + pkin(3);
t36 = qJ(4) + rSges(5,3);
t30 = qJ(1) + qJ(2);
t27 = sin(t30);
t35 = pkin(2) * t27;
t31 = sin(qJ(1));
t34 = t31 * pkin(1);
t33 = Icges(3,3) + Icges(4,3) + Icges(5,2);
t28 = cos(t30);
t12 = t28 * rSges(3,1) - t27 * rSges(3,2);
t26 = pkin(6) + t30;
t23 = sin(t26);
t24 = cos(t26);
t25 = pkin(2) * t28;
t8 = t24 * rSges(4,1) - t23 * rSges(4,2) + t25;
t4 = t36 * t23 + t37 * t24 + t25;
t11 = -t27 * rSges(3,1) - t28 * rSges(3,2);
t7 = -t23 * rSges(4,1) - t24 * rSges(4,2) - t35;
t3 = -t37 * t23 + t36 * t24 - t35;
t32 = cos(qJ(1));
t29 = t32 * pkin(1);
t14 = t32 * rSges(2,1) - t31 * rSges(2,2);
t13 = -t31 * rSges(2,1) - t32 * rSges(2,2);
t10 = t12 + t29;
t9 = t11 - t34;
t6 = t29 + t8;
t5 = t7 - t34;
t2 = t29 + t4;
t1 = t3 - t34;
t15 = [Icges(2,3) + m(2) * (t13 ^ 2 + t14 ^ 2) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t33; m(3) * (t12 * t10 + t11 * t9) + m(4) * (t7 * t5 + t8 * t6) + m(5) * (t3 * t1 + t4 * t2) + t33; m(3) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t33; 0; 0; m(4) + m(5); m(5) * (t23 * t1 - t24 * t2); m(5) * (t23 * t3 - t24 * t4); 0; m(5) * (t23 ^ 2 + t24 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;
