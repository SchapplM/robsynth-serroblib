% Calculate joint inertia matrix for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:25
% EndTime: 2019-03-08 18:27:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (202->47), mult. (192->60), div. (0->0), fcn. (160->6), ass. (0->24)
t23 = sin(qJ(1));
t29 = t23 * pkin(1);
t28 = cos(qJ(4));
t27 = sin(qJ(4));
t22 = qJ(1) + pkin(6);
t19 = sin(t22);
t20 = cos(t22);
t24 = cos(qJ(1));
t21 = t24 * pkin(1);
t26 = t20 * pkin(2) + t19 * qJ(3) + t21;
t25 = t20 * qJ(3) - t29;
t10 = t19 * t28 - t20 * t27;
t9 = -t19 * t27 - t20 * t28;
t4 = t9 * rSges(5,1) - t10 * rSges(5,2);
t3 = t10 * rSges(5,1) + t9 * rSges(5,2);
t13 = t24 * rSges(2,1) - t23 * rSges(2,2);
t12 = -t23 * rSges(2,1) - t24 * rSges(2,2);
t8 = t20 * rSges(3,1) - t19 * rSges(3,2) + t21;
t7 = -t19 * rSges(3,1) - t20 * rSges(3,2) - t29;
t6 = t20 * rSges(4,1) + t19 * rSges(4,3) + t26;
t5 = t20 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t19 + t25;
t2 = t20 * pkin(3) + t26 - t4;
t1 = (-pkin(2) - pkin(3)) * t19 - t3 + t25;
t11 = [Icges(2,3) + Icges(3,3) + Icges(4,2) + Icges(5,3) + m(2) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; m(3) + m(4) + m(5); m(4) * (t19 * t5 - t20 * t6) + m(5) * (t19 * t1 - t20 * t2); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t19 ^ 2 + t20 ^ 2); -Icges(5,3) + m(5) * (t3 * t1 + t4 * t2); 0; m(5) * (t3 * t19 - t4 * t20); m(5) * (t3 ^ 2 + t4 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
