% Calculate joint inertia matrix for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (192->59), mult. (344->77), div. (0->0), fcn. (312->4), ass. (0->28)
t32 = cos(qJ(3));
t23 = sin(qJ(3));
t24 = sin(qJ(1));
t31 = t24 * t23;
t25 = cos(qJ(1));
t30 = t25 * t23;
t29 = t25 * pkin(1) + t24 * qJ(2);
t28 = Icges(4,3) + Icges(5,3);
t11 = -t25 * t32 - t31;
t12 = t24 * t32 - t30;
t7 = t12 * rSges(4,1) + t11 * rSges(4,2);
t8 = t11 * rSges(4,1) - t12 * rSges(4,2);
t27 = t12 * rSges(5,1) + t11 * rSges(5,2) - pkin(3) * t30;
t18 = t32 * pkin(3) + pkin(2);
t26 = t11 * rSges(5,1) - t12 * rSges(5,2) - pkin(3) * t31 - t25 * t18;
t21 = t25 * pkin(2);
t20 = t25 * qJ(2);
t14 = t25 * rSges(2,1) - t24 * rSges(2,2);
t13 = -t24 * rSges(2,1) - t25 * rSges(2,2);
t10 = t25 * rSges(3,1) + t24 * rSges(3,3) + t29;
t9 = t25 * rSges(3,3) + t20 + (-rSges(3,1) - pkin(1)) * t24;
t6 = t21 - t8 + t29;
t5 = t20 + (-pkin(1) - pkin(2)) * t24 - t7;
t4 = t21 + t26;
t3 = (-pkin(2) + t18) * t24 + t27;
t2 = -t26 + t29;
t1 = t20 + (-pkin(1) - t18) * t24 - t27;
t15 = [Icges(2,3) + Icges(3,2) + m(2) * (t13 ^ 2 + t14 ^ 2) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t28; m(3) * (-t25 * t10 + t24 * t9) + m(4) * (t24 * t5 - t25 * t6) + m(5) * (t24 * t1 - t25 * t2); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t24 ^ 2 + t25 ^ 2); m(4) * (t7 * t5 + t8 * t6) + m(5) * (t3 * t1 + t4 * t2) - t28; m(4) * (t7 * t24 - t8 * t25) + m(5) * (t3 * t24 - t4 * t25); m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t28; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;
