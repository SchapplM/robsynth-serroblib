% Calculate joint inertia matrix for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:19
% EndTime: 2018-11-14 14:01:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (58->32), mult. (112->40), div. (0->0), fcn. (61->2), ass. (0->13)
t15 = -m(4) - m(5);
t14 = rSges(5,1) + pkin(3);
t11 = sin(qJ(2));
t12 = cos(qJ(2));
t13 = t12 * pkin(2) + t11 * qJ(3);
t9 = t12 * qJ(3);
t6 = rSges(3,1) * t12 - t11 * rSges(3,2);
t5 = -t11 * rSges(3,1) - rSges(3,2) * t12;
t4 = rSges(4,1) * t12 + t11 * rSges(4,3) + t13;
t3 = rSges(4,3) * t12 + t9 + (-rSges(4,1) - pkin(2)) * t11;
t2 = t11 * rSges(5,2) + t12 * t14 + t13;
t1 = rSges(5,2) * t12 + t9 + (-pkin(2) - t14) * t11;
t7 = [m(2) + m(3) - t15; m(3) * t6 + m(4) * t4 + m(5) * t2; Icges(3,3) + Icges(4,2) + Icges(5,3) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); t15 * t12; m(4) * (t11 * t3 - t12 * t4) + m(5) * (t11 * t1 - t12 * t2); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t11 ^ 2 + t12 ^ 2); 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
