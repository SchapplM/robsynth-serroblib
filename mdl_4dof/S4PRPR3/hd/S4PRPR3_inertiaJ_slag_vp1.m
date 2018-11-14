% Calculate joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:15
% EndTime: 2018-11-14 14:11:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (98->31), mult. (94->35), div. (0->0), fcn. (45->6), ass. (0->20)
t20 = m(4) + m(5);
t17 = sin(qJ(2));
t19 = t17 * pkin(2);
t16 = qJ(2) + pkin(6);
t14 = qJ(4) + t16;
t10 = sin(t14);
t11 = cos(t14);
t6 = t11 * rSges(5,1) - t10 * rSges(5,2);
t5 = -t10 * rSges(5,1) - t11 * rSges(5,2);
t18 = cos(qJ(2));
t15 = t18 * pkin(2);
t13 = cos(t16);
t12 = sin(t16);
t8 = t18 * rSges(3,1) - t17 * rSges(3,2);
t7 = -t17 * rSges(3,1) - t18 * rSges(3,2);
t4 = t13 * rSges(4,1) - t12 * rSges(4,2) + t15;
t3 = -t12 * rSges(4,1) - t13 * rSges(4,2) - t19;
t2 = pkin(3) * t13 + t15 + t6;
t1 = -pkin(3) * t12 - t19 + t5;
t9 = [m(2) + m(3) + t20; m(3) * t8 + m(4) * t4 + m(5) * t2; Icges(3,3) + Icges(4,3) + Icges(5,3) + m(3) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; 0; t20; m(5) * t6; Icges(5,3) + m(5) * (t5 * t1 + t6 * t2); 0; m(5) * (t5 ^ 2 + t6 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1) t9(2) t9(4) t9(7); t9(2) t9(3) t9(5) t9(8); t9(4) t9(5) t9(6) t9(9); t9(7) t9(8) t9(9) t9(10);];
Mq  = res;
