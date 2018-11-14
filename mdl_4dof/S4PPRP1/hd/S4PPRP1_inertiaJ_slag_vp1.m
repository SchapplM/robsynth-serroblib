% Calculate joint inertia matrix for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:04
% EndTime: 2018-11-14 13:39:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (79->25), mult. (168->40), div. (0->0), fcn. (164->4), ass. (0->13)
t17 = rSges(5,1) + pkin(3);
t16 = cos(qJ(3));
t15 = sin(qJ(3));
t14 = -rSges(5,3) - qJ(4);
t13 = cos(pkin(5));
t12 = sin(pkin(5));
t6 = -t12 * t16 + t13 * t15;
t5 = -t12 * t15 - t13 * t16;
t4 = t5 * rSges(4,1) + t6 * rSges(4,2);
t3 = -t6 * rSges(4,1) + t5 * rSges(4,2);
t2 = t14 * t6 + t17 * t5;
t1 = t14 * t5 - t17 * t6;
t7 = [m(2) + m(3) + m(4) + m(5); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t12 ^ 2 + t13 ^ 2); 0; m(4) * (t3 * t12 - t4 * t13) + m(5) * (t1 * t12 - t2 * t13); Icges(4,3) + Icges(5,2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; m(5) * (t6 * t12 + t5 * t13); m(5) * (t6 * t1 - t5 * t2); m(5) * (t5 ^ 2 + t6 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
