% Calculate joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:06
% DurationCPUTime: 0.68s
% Computational Cost: add. (694->108), mult. (1388->148), div. (0->0), fcn. (1568->6), ass. (0->51)
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t79 = sin(pkin(7));
t80 = cos(pkin(7));
t35 = -t59 * t80 + t61 * t79;
t32 = t35 ^ 2;
t34 = -t59 * t79 - t61 * t80;
t33 = t34 ^ 2;
t86 = t32 + t33;
t107 = Icges(5,5) + Icges(6,5);
t106 = Icges(5,6) + Icges(6,6);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t104 = t106 * t58 - t107 * t60;
t101 = Icges(5,3) + Icges(6,3);
t88 = rSges(6,3) + qJ(5) + pkin(6);
t98 = t101 * t34 - t104 * t35;
t97 = t101 * t35 + t104 * t34;
t45 = -t58 * rSges(5,1) - t60 * rSges(5,2);
t94 = m(5) * t45;
t93 = rSges(5,2) * t58;
t92 = rSges(6,2) * t58;
t91 = t34 * rSges(5,3);
t90 = t34 * t58;
t89 = t34 * t60;
t87 = -rSges(5,1) * t89 + t35 * rSges(5,3);
t85 = t61 * pkin(1) + t59 * qJ(2);
t78 = t61 * pkin(2) + t85;
t77 = t60 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t58;
t52 = t60 * pkin(4) + pkin(3);
t76 = -rSges(6,1) * t89 - t34 * t52 + t88 * t35;
t54 = t61 * qJ(2);
t74 = t54 + (-pkin(1) - pkin(2)) * t59;
t73 = -rSges(5,1) * t60 + t93;
t62 = rSges(6,1) * t60 + t52 - t92;
t47 = t61 * rSges(2,1) - t59 * rSges(2,2);
t46 = -t59 * rSges(2,1) - t61 * rSges(2,2);
t31 = t34 * pkin(6);
t25 = t61 * rSges(3,1) + t59 * rSges(3,3) + t85;
t24 = t61 * rSges(3,3) + t54 + (-rSges(3,1) - pkin(1)) * t59;
t22 = t77 * t34;
t21 = t77 * t35;
t20 = -t34 * rSges(4,1) - t35 * rSges(4,2) + t78;
t19 = t35 * rSges(4,1) - t34 * rSges(4,2) + t74;
t6 = t35 * pkin(6) + (-pkin(3) + t93) * t34 + t78 + t87;
t5 = t91 + t31 + (pkin(3) - t73) * t35 + t74;
t4 = rSges(6,2) * t90 + t76 + t78;
t3 = t34 * t88 + t62 * t35 + t74;
t2 = t34 * (rSges(5,2) * t90 + t87) + (t35 * t73 - t91) * t35;
t1 = ((pkin(3) + t92) * t34 + t76) * t34 + (t31 + (pkin(3) - t62) * t35 + (-pkin(6) - t88) * t34) * t35;
t7 = [Icges(3,2) + Icges(2,3) + Icges(4,3) + m(2) * (t46 ^ 2 + t47 ^ 2) + m(3) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (Icges(5,2) + Icges(6,2)) * t60 ^ 2 + (0.2e1 * (Icges(6,4) + Icges(5,4)) * t60 + (Icges(5,1) + Icges(6,1)) * t58) * t58; m(3) * (t59 * t24 - t61 * t25) + m(4) * (t59 * t19 - t61 * t20) + m(5) * (t59 * t5 - t61 * t6) + m(6) * (t59 * t3 - t61 * t4); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t59 ^ 2 + t61 ^ 2); 0; 0; m(4) + m(5) + m(6); m(6) * (t21 * t4 + t22 * t3) + (-t34 * t5 - t35 * t6) * t94 - t86 * (t106 * t60 + t107 * t58); m(6) * (-t21 * t61 + t22 * t59) + (-t34 * t59 + t35 * t61) * t94; -m(5) * t2 - m(6) * t1; m(5) * (t45 ^ 2 * t86 + t2 ^ 2) + m(6) * (t1 ^ 2 + t21 ^ 2 + t22 ^ 2) + t97 * t35 * t32 + (t98 * t33 + (t97 * t34 + t98 * t35) * t35) * t34; m(6) * (t35 * t3 - t34 * t4); m(6) * (t34 * t61 + t35 * t59); 0; m(6) * (-t34 * t21 + t35 * t22); m(6) * t86;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
