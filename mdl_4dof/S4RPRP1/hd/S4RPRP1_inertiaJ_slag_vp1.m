% Calculate joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (250->44), mult. (160->53), div. (0->0), fcn. (100->6), ass. (0->28)
t35 = rSges(5,1) + pkin(3);
t34 = qJ(4) + rSges(5,3);
t28 = sin(qJ(1));
t33 = t28 * pkin(1);
t27 = qJ(1) + pkin(6);
t24 = cos(t27);
t29 = cos(qJ(1));
t26 = t29 * pkin(1);
t32 = pkin(2) * t24 + t26;
t31 = Icges(4,3) + Icges(5,2);
t25 = qJ(3) + t27;
t21 = sin(t25);
t22 = cos(t25);
t10 = t22 * rSges(4,1) - t21 * rSges(4,2);
t4 = t34 * t21 + t35 * t22;
t23 = sin(t27);
t30 = -pkin(2) * t23 - t33;
t9 = -t21 * rSges(4,1) - t22 * rSges(4,2);
t3 = -t35 * t21 + t34 * t22;
t12 = t29 * rSges(2,1) - t28 * rSges(2,2);
t11 = -t28 * rSges(2,1) - t29 * rSges(2,2);
t8 = t24 * rSges(3,1) - t23 * rSges(3,2) + t26;
t7 = -t23 * rSges(3,1) - t24 * rSges(3,2) - t33;
t6 = t10 + t32;
t5 = t30 + t9;
t2 = t4 + t32;
t1 = t3 + t30;
t13 = [Icges(2,3) + Icges(3,3) + m(2) * (t11 ^ 2 + t12 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t31; 0; m(3) + m(4) + m(5); m(4) * (t10 * t6 + t9 * t5) + m(5) * (t3 * t1 + t4 * t2) + t31; 0; m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t31; m(5) * (t21 * t1 - t22 * t2); 0; m(5) * (t21 * t3 - t22 * t4); m(5) * (t21 ^ 2 + t22 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1) t13(2) t13(4) t13(7); t13(2) t13(3) t13(5) t13(8); t13(4) t13(5) t13(6) t13(9); t13(7) t13(8) t13(9) t13(10);];
Mq  = res;
