% Calculate time derivative of joint inertia matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.29s
% Computational Cost: add. (1317->90), mult. (1164->125), div. (0->0), fcn. (892->6), ass. (0->58)
t57 = qJD(1) + qJD(2);
t79 = -qJD(4) + t57;
t58 = qJ(1) + qJ(2);
t54 = sin(t58);
t55 = cos(t58);
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t78 = -t54 * t61 + t55 * t59;
t77 = -pkin(2) - pkin(3);
t60 = sin(qJ(1));
t76 = t60 * pkin(1);
t75 = -rSges(4,1) - pkin(2);
t74 = t54 * rSges(3,2);
t72 = t55 * t57;
t63 = t54 * t59 + t55 * t61;
t17 = t79 * t63;
t18 = t79 * t78;
t6 = -t18 * rSges(5,1) - t17 * rSges(5,2);
t20 = -rSges(5,1) * t63 + rSges(5,2) * t78;
t47 = t55 * qJ(3);
t70 = qJD(3) * t54 + t57 * t47;
t46 = t54 * qJ(3);
t69 = t55 * pkin(2) + t46;
t68 = pkin(1) * qJD(1);
t67 = t60 * t68;
t62 = cos(qJ(1));
t66 = t62 * t68;
t65 = t77 * t54;
t64 = t75 * t54;
t36 = t55 * rSges(3,1) - t74;
t24 = t55 * rSges(4,1) + t54 * rSges(4,3) + t69;
t28 = -rSges(3,1) * t72 + t57 * t74;
t10 = t55 * pkin(3) - t20 + t69;
t35 = -t54 * rSges(3,1) - t55 * rSges(3,2);
t5 = t17 * rSges(5,1) - t18 * rSges(5,2);
t19 = -rSges(5,1) * t78 - rSges(5,2) * t63;
t23 = t55 * rSges(4,3) + t47 + t64;
t27 = t35 * t57;
t13 = rSges(4,3) * t72 + t57 * t64 + t70;
t3 = t57 * t65 - t6 + t70;
t9 = t47 + t65 - t19;
t45 = qJD(3) * t55;
t14 = t45 + (t75 * t55 + (-rSges(4,3) - qJ(3)) * t54) * t57;
t4 = t45 + (t77 * t55 - t46) * t57 - t5;
t56 = t62 * pkin(1);
t32 = t36 + t56;
t31 = t35 - t76;
t26 = t28 - t66;
t25 = t27 - t67;
t22 = t56 + t24;
t21 = t23 - t76;
t12 = t14 - t66;
t11 = t13 - t67;
t8 = t56 + t10;
t7 = t9 - t76;
t2 = t4 - t66;
t1 = t3 - t67;
t15 = [0.2e1 * m(3) * (t32 * t25 + t31 * t26) + 0.2e1 * m(4) * (t22 * t11 + t21 * t12) + 0.2e1 * m(5) * (t8 * t1 + t7 * t2); m(3) * (t36 * t25 + t35 * t26 + t27 * t32 + t28 * t31) + m(4) * (t24 * t11 + t23 * t12 + t13 * t22 + t14 * t21) + m(5) * (t10 * t1 + t9 * t2 + t3 * t8 + t4 * t7); 0.2e1 * m(3) * (t36 * t27 + t35 * t28) + 0.2e1 * m(4) * (t24 * t13 + t23 * t14) + 0.2e1 * m(5) * (t10 * t3 + t9 * t4); m(4) * ((t21 * t57 - t11) * t55 + (t22 * t57 + t12) * t54) + m(5) * ((t57 * t7 - t1) * t55 + (t57 * t8 + t2) * t54); m(4) * ((t23 * t57 - t13) * t55 + (t24 * t57 + t14) * t54) + m(5) * ((t57 * t9 - t3) * t55 + (t10 * t57 + t4) * t54); 0; m(5) * (t20 * t1 + t19 * t2 + t5 * t7 + t6 * t8); m(5) * (t6 * t10 + t19 * t4 + t20 * t3 + t5 * t9); m(5) * ((t19 * t57 - t6) * t55 + (t20 * t57 + t5) * t54); 0.2e1 * m(5) * (t19 * t5 + t20 * t6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1), t15(2), t15(4), t15(7); t15(2), t15(3), t15(5), t15(8); t15(4), t15(5), t15(6), t15(9); t15(7), t15(8), t15(9), t15(10);];
Mq  = res;
