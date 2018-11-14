% Calculate time derivative of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:26
% DurationCPUTime: 0.34s
% Computational Cost: add. (1161->109), mult. (1694->153), div. (0->0), fcn. (1354->9), ass. (0->62)
t38 = sin(pkin(6));
t41 = sin(qJ(1));
t42 = cos(qJ(1));
t55 = pkin(4) + pkin(6);
t56 = pkin(4) - pkin(6);
t46 = cos(t56) / 0.2e1 + cos(t55) / 0.2e1;
t43 = t42 * t46;
t24 = t41 * t38 - t43;
t39 = sin(pkin(4));
t64 = qJ(2) * t39;
t36 = t42 * t64;
t75 = -t24 * qJ(3) + t36;
t74 = 0.2e1 * t39;
t73 = m(4) / 0.2e1;
t72 = m(5) / 0.2e1;
t71 = t41 * pkin(1);
t37 = t42 * pkin(1);
t70 = rSges(5,1) + pkin(3);
t69 = rSges(4,2) - pkin(2);
t68 = t39 * t41;
t67 = t39 * t42;
t58 = qJD(2) * t39;
t66 = qJD(1) * t36 + t41 * t58;
t65 = t41 * t64 + t37;
t63 = rSges(5,2) + qJ(3);
t62 = rSges(4,3) + qJ(3);
t61 = rSges(5,3) + qJ(4);
t60 = qJD(1) * t41;
t59 = qJD(1) * t42;
t57 = t73 + t72;
t54 = -pkin(2) - t61;
t40 = cos(pkin(6));
t45 = sin(t55) / 0.2e1 - sin(t56) / 0.2e1;
t44 = t41 * t45;
t27 = t42 * t40 - t44;
t53 = t27 * pkin(2) + t65;
t26 = t42 * t38 + t41 * t46;
t19 = t26 * qJD(1);
t34 = t42 * t58;
t52 = -t19 * qJ(3) - t24 * qJD(3) + t34;
t51 = t39 * t70;
t25 = t41 * t40 + t42 * t45;
t18 = t25 * qJD(1);
t50 = -t18 * pkin(2) + t26 * qJD(3) + t66;
t49 = rSges(4,1) * t67 - t71;
t48 = rSges(3,3) * t67 - t71;
t47 = t42 * t51 - t71;
t20 = -qJD(1) * t44 + t40 * t59;
t17 = -qJD(1) * t43 + t38 * t60;
t14 = rSges(3,1) * t27 - rSges(3,2) * t26 + rSges(3,3) * t68 + t65;
t13 = -t25 * rSges(3,1) + t24 * rSges(3,2) + t36 + t48;
t12 = -t20 * rSges(3,1) + t19 * rSges(3,2) + t34 + (-t37 + (-rSges(3,3) - qJ(2)) * t68) * qJD(1);
t11 = -t18 * rSges(3,1) + t17 * rSges(3,2) + t48 * qJD(1) + t66;
t9 = rSges(4,1) * t68 - t27 * rSges(4,2) + t62 * t26 + t53;
t8 = -t24 * rSges(4,3) + t69 * t25 + t49 + t75;
t6 = t63 * t26 + t61 * t27 + t41 * t51 + t53;
t5 = -t24 * rSges(5,2) + t54 * t25 + t47 + t75;
t4 = -t19 * rSges(4,3) + t69 * t20 + (-t37 + (-rSges(4,1) - qJ(2)) * t68) * qJD(1) + t52;
t3 = t18 * rSges(4,2) + t49 * qJD(1) - t62 * t17 + t50;
t2 = -t19 * rSges(5,2) - t25 * qJD(4) + t54 * t20 + (-t37 + (-qJ(2) - t70) * t68) * qJD(1) + t52;
t1 = t47 * qJD(1) + t27 * qJD(4) - t63 * t17 - t61 * t18 + t50;
t7 = [0.2e1 * m(3) * (t11 * t14 + t12 * t13) + 0.2e1 * m(4) * (t3 * t9 + t4 * t8) + 0.2e1 * m(5) * (t1 * t6 + t2 * t5); (m(3) * (-t11 * t42 + t12 * t41 + t13 * t59 + t14 * t60) / 0.2e1 + (-t3 * t42 + t4 * t41 + t8 * t59 + t9 * t60) * t73 + (-t1 * t42 + t2 * t41 + t5 * t59 + t6 * t60) * t72) * t74; 0; m(4) * (-t17 * t8 + t19 * t9 + t24 * t3 + t26 * t4) + m(5) * (t1 * t24 - t17 * t5 + t19 * t6 + t2 * t26); t57 * (-t17 * t41 - t19 * t42 + (t24 * t41 + t26 * t42) * qJD(1)) * t74; 0.4e1 * t57 * (-t17 * t26 + t19 * t24); m(5) * (t1 * t25 - t18 * t5 + t2 * t27 + t20 * t6); m(5) * (-t18 * t41 - t20 * t42 + (t25 * t41 + t27 * t42) * qJD(1)) * t39; m(5) * (-t17 * t27 - t18 * t26 + t19 * t25 + t20 * t24); 0.2e1 * m(5) * (-t18 * t27 + t20 * t25);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
