% Calculate time derivative of joint inertia matrix for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S2RR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_inertiaDJ_slag_vp1: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_inertiaDJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR2_inertiaDJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR2_inertiaDJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:45
% EndTime: 2018-11-16 16:48:46
% DurationCPUTime: 0.65s
% Computational Cost: add. (370->94), mult. (1022->161), div. (0->0), fcn. (786->4), ass. (0->65)
t38 = cos(qJ(2));
t36 = sin(qJ(2));
t75 = rSges(3,1) * t36;
t30 = rSges(3,2) * t38 + t75;
t39 = cos(qJ(1));
t56 = qJD(2) * t39;
t84 = t30 * t56;
t37 = sin(qJ(1));
t65 = Icges(3,4) * t38;
t46 = -Icges(3,2) * t36 + t65;
t17 = Icges(3,6) * t39 - t46 * t37;
t66 = Icges(3,4) * t36;
t48 = Icges(3,1) * t38 - t66;
t19 = Icges(3,5) * t39 - t48 * t37;
t50 = t17 * t36 - t19 * t38;
t83 = t50 * t39;
t73 = rSges(3,2) * t36;
t74 = rSges(3,1) * t38;
t52 = -t73 + t74;
t82 = t52 * t39;
t72 = rSges(3,2) * t37;
t81 = t39 * rSges(3,3) + t36 * t72;
t44 = Icges(3,5) * t38 - Icges(3,6) * t36;
t16 = Icges(3,3) * t37 + t44 * t39;
t18 = Icges(3,6) * t37 + t46 * t39;
t20 = Icges(3,5) * t37 + t48 * t39;
t80 = 2 * m(3);
t79 = t37 ^ 2;
t78 = t39 ^ 2;
t77 = m(3) * t30;
t76 = -rSges(3,3) - pkin(1);
t71 = t17 * t38;
t70 = t18 * t38;
t69 = t19 * t36;
t68 = t20 * t36;
t67 = t37 * rSges(3,3);
t15 = Icges(3,3) * t39 - t44 * t37;
t61 = qJD(1) * t15;
t60 = qJD(1) * t39;
t59 = qJD(2) * t36;
t58 = qJD(2) * t37;
t57 = qJD(2) * t38;
t55 = -t57 * t72 - t58 * t75 - t60 * t73;
t54 = t76 * t37;
t49 = t18 * t36 - t20 * t38;
t47 = Icges(3,1) * t36 + t65;
t45 = Icges(3,2) * t38 + t66;
t43 = Icges(3,5) * t36 + Icges(3,6) * t38;
t21 = -t37 * t74 + t81;
t42 = t49 * t37;
t41 = qJD(2) * t47;
t40 = qJD(2) * t45;
t26 = t52 * qJD(2);
t22 = t67 + t82;
t14 = pkin(1) * t39 + t21;
t13 = t54 - t82;
t8 = -qJD(1) * t16 + t43 * t58;
t7 = -t43 * t56 + t61;
t6 = t84 + (t52 * t37 + t76 * t39) * qJD(1);
t5 = (-t39 * t74 + t54) * qJD(1) - t55;
t4 = t37 * t16 - t49 * t39;
t3 = t37 * t15 - t83;
t2 = t16 * t39 + t42;
t1 = t15 * t39 + t50 * t37;
t9 = [(t13 * t6 + t14 * t5) * t80 + t38 * t41 + t48 * t59 - t36 * t40 + t46 * t57; m(3) * ((t37 * t6 - t39 * t5) * t30 + (t13 * t37 - t14 * t39) * t26) + (-t49 * qJD(2) + t36 * (t19 * qJD(1) - t47 * t56) + t38 * (t17 * qJD(1) - t45 * t56)) * t37 / 0.2e1 + (-t50 * qJD(2) + t38 * (-t18 * qJD(1) + t37 * t40) + t36 * (-t20 * qJD(1) + t37 * t41)) * t39 / 0.2e1 + (t79 / 0.2e1 + t78 / 0.2e1) * t44 * qJD(2) + ((t14 * t77 - t71 / 0.2e1 - t69 / 0.2e1) * t37 + (t13 * t77 + t70 / 0.2e1 + t68 / 0.2e1) * t39) * qJD(1); ((-t37 * t21 + t22 * t39) * (((-t22 + t67) * qJD(1) + t55) * t37 + (-t84 + (-t21 + t81) * qJD(1)) * t39) + (t78 + t79) * t30 * t26) * t80 + (t3 * t39 + t4 * t37) * t60 + t37 * ((t37 * t7 + (-t3 + t42) * qJD(1)) * t37 + (t4 * qJD(1) + (-t17 * t57 - t19 * t59 + t61) * t39 + (t8 + (-t68 - t70) * qJD(2) + t50 * qJD(1)) * t37) * t39) - qJD(1) * t37 * (t1 * t39 + t2 * t37) + t39 * ((t39 * t8 + (t2 + t83) * qJD(1)) * t39 + (-t1 * qJD(1) + (t18 * t57 + t20 * t59) * t37 + (t7 + (t69 + t71) * qJD(2) + (-t15 + t49) * qJD(1)) * t39) * t37);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t9(1) t9(2); t9(2) t9(3);];
Mq  = res;
