% Calculate joint inertia matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:45
% DurationCPUTime: 0.47s
% Computational Cost: add. (543->72), mult. (562->97), div. (0->0), fcn. (484->4), ass. (0->39)
t43 = pkin(6) + qJ(2);
t41 = sin(t43);
t39 = t41 ^ 2;
t42 = cos(t43);
t40 = t42 ^ 2;
t66 = t39 + t40;
t85 = Icges(4,5) + Icges(5,5);
t84 = Icges(4,6) + Icges(5,6);
t45 = sin(qJ(3));
t46 = cos(qJ(3));
t82 = -t84 * t45 + t85 * t46;
t79 = Icges(4,3) + Icges(5,3);
t76 = -t82 * t41 + t79 * t42;
t75 = t79 * t41 + t82 * t42;
t72 = rSges(4,2) * t45;
t71 = rSges(5,2) * t45;
t70 = t42 * rSges(4,3);
t69 = t42 * t45;
t68 = t42 * t46;
t67 = rSges(4,1) * t68 + t41 * rSges(4,3);
t38 = t46 * pkin(3) + pkin(2);
t61 = rSges(5,1) * t68 + t41 * rSges(5,3) + t42 * t38;
t60 = -t46 * rSges(5,2) + (-rSges(5,1) - pkin(3)) * t45;
t58 = rSges(4,1) * t46 - t72;
t47 = -rSges(5,1) * t46 - t38 + t71;
t44 = -qJ(4) - pkin(5);
t37 = t42 * pkin(5);
t32 = t45 * rSges(4,1) + t46 * rSges(4,2);
t22 = t42 * rSges(3,1) - t41 * rSges(3,2);
t21 = -t41 * rSges(3,1) - t42 * rSges(3,2);
t20 = t60 * t42;
t19 = t60 * t41;
t6 = t41 * pkin(5) + (pkin(2) - t72) * t42 + t67;
t5 = t70 + t37 + (-pkin(2) - t58) * t41;
t4 = -rSges(5,2) * t69 - t41 * t44 + t61;
t3 = (rSges(5,3) - t44) * t42 + t47 * t41;
t2 = t42 * (-rSges(4,2) * t69 + t67) + (t58 * t41 - t70) * t41;
t1 = ((-pkin(2) - t71) * t42 + t61) * t42 + (t37 + (-pkin(2) - t47) * t41 + (-pkin(5) - rSges(5,3)) * t42) * t41;
t7 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,2) + Icges(5,2)) * t46 ^ 2 + (0.2e1 * (Icges(5,4) + Icges(4,4)) * t46 + (Icges(4,1) + Icges(5,1)) * t45) * t45; m(4) * t2 + m(5) * t1; m(5) * (t19 * t4 + t20 * t3) + m(4) * (-t41 * t6 - t42 * t5) * t32 + t66 * (t85 * t45 + t84 * t46); m(4) * (t66 * t32 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t75 * t41 * t39 + (t76 * t40 + (t76 * t41 + t75 * t42) * t41) * t42; 0; m(5) * (t41 * t3 - t42 * t4); m(5) * (-t42 * t19 + t41 * t20); m(5) * t66;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
