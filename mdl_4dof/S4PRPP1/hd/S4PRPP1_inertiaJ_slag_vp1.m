% Calculate joint inertia matrix for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:54
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (119->32), mult. (104->37), div. (0->0), fcn. (62->2), ass. (0->15)
t14 = pkin(5) + qJ(2);
t12 = sin(t14);
t13 = cos(t14);
t16 = t13 * pkin(2) + t12 * qJ(3);
t15 = rSges(5,3) + qJ(4);
t10 = t13 * qJ(3);
t8 = t12 ^ 2 + t13 ^ 2;
t7 = rSges(3,1) * t13 - rSges(3,2) * t12;
t6 = -rSges(3,1) * t12 - rSges(3,2) * t13;
t5 = m(5) * t8;
t4 = -rSges(4,2) * t13 + rSges(4,3) * t12 + t16;
t3 = rSges(4,3) * t13 + t10 + (rSges(4,2) - pkin(2)) * t12;
t2 = rSges(5,2) * t12 + t15 * t13 + t16;
t1 = rSges(5,2) * t13 + t10 + (-pkin(2) - t15) * t12;
t9 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + Icges(4,1) + Icges(5,1) + m(3) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; m(4) * (t12 * t3 - t13 * t4) + m(5) * (t12 * t1 - t13 * t2); m(4) * t8 + t5; 0; m(5) * (t1 * t13 + t12 * t2); 0; t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1) t9(2) t9(4) t9(7); t9(2) t9(3) t9(5) t9(8); t9(4) t9(5) t9(6) t9(9); t9(7) t9(8) t9(9) t9(10);];
Mq  = res;
