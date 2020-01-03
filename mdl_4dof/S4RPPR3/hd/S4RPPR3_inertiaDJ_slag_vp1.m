% Calculate time derivative of joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.81s
% Computational Cost: add. (1448->121), mult. (1322->181), div. (0->0), fcn. (986->8), ass. (0->78)
t48 = pkin(7) + qJ(4);
t43 = sin(t48);
t106 = qJD(1) * t43;
t45 = cos(t48);
t105 = qJD(1) * t45;
t49 = qJ(1) + pkin(6);
t44 = sin(t49);
t35 = rSges(5,1) * t43 + rSges(5,2) * t45;
t60 = t35 * qJD(4);
t104 = t44 * t60;
t46 = cos(t49);
t84 = Icges(5,4) * t45;
t65 = -Icges(5,2) * t43 + t84;
t23 = Icges(5,6) * t44 + t65 * t46;
t85 = Icges(5,4) * t43;
t67 = Icges(5,1) * t45 - t85;
t25 = Icges(5,5) * t44 + t67 * t46;
t68 = t23 * t43 - t25 * t45;
t103 = t68 * t44;
t22 = -Icges(5,6) * t46 + t65 * t44;
t24 = -Icges(5,5) * t46 + t67 * t44;
t69 = t22 * t43 - t24 * t45;
t102 = t69 * t46;
t51 = cos(pkin(7));
t42 = pkin(3) * t51 + pkin(2);
t94 = rSges(5,2) * t43;
t95 = rSges(5,1) * t45;
t72 = -t94 + t95;
t61 = -t42 - t72;
t52 = -pkin(5) - qJ(3);
t87 = rSges(5,3) - t52;
t96 = sin(qJ(1)) * pkin(1);
t16 = t61 * t44 + t87 * t46 - t96;
t47 = cos(qJ(1)) * pkin(1);
t41 = t44 * rSges(5,3);
t86 = t46 * t95 + t41;
t17 = -t44 * t52 + t47 + (t42 - t94) * t46 + t86;
t101 = t16 * t46 + t17 * t44;
t63 = Icges(5,5) * t45 - Icges(5,6) * t43;
t20 = -Icges(5,3) * t46 + t63 * t44;
t62 = rSges(4,1) * t51 - rSges(4,2) * sin(pkin(7)) + pkin(2);
t80 = rSges(4,3) + qJ(3);
t19 = t80 * t44 + t62 * t46 + t47;
t100 = 2 * m(5);
t99 = t44 ^ 2;
t98 = t46 ^ 2;
t97 = m(5) * t35;
t91 = t43 * t24;
t90 = t43 * t25;
t89 = t45 * t22;
t88 = t45 * t23;
t21 = Icges(5,3) * t44 + t63 * t46;
t79 = qJD(1) * t21;
t78 = qJD(1) * t44;
t77 = qJD(1) * t46;
t76 = qJD(4) * t43;
t75 = qJD(4) * t45;
t74 = t46 * t94;
t57 = qJD(4) * (-Icges(5,5) * t43 - Icges(5,6) * t45);
t55 = rSges(5,3) * t77 - t46 * t60 + t78 * t94;
t18 = -t62 * t44 + t80 * t46 - t96;
t40 = qJD(3) * t46;
t39 = qJD(3) * t44;
t31 = t72 * qJD(4);
t27 = -t74 + t86;
t26 = -rSges(5,3) * t46 + t72 * t44;
t15 = -qJD(1) * t19 + t40;
t14 = qJD(1) * t18 + t39;
t9 = t44 * t57 + t79;
t8 = -qJD(1) * t20 + t46 * t57;
t7 = t40 + t104 + (-t87 * t44 + t61 * t46 - t47) * qJD(1);
t6 = t39 + (-t96 - t46 * t52 + (-t42 - t95) * t44) * qJD(1) + t55;
t5 = t21 * t44 - t68 * t46;
t4 = t20 * t44 - t102;
t3 = -t21 * t46 - t103;
t2 = -t20 * t46 - t69 * t44;
t1 = (qJD(1) * t26 + t55) * t46 + (-t104 + (-t27 - t74 + t41) * qJD(1)) * t44;
t10 = [0.2e1 * m(4) * (t14 * t19 + t15 * t18) + (t16 * t7 + t17 * t6) * t100 + (-Icges(5,2) * t45 + t67 - t85) * t76 + (Icges(5,1) * t43 + t65 + t84) * t75; 0; 0; m(4) * (-t14 * t46 + t15 * t44 + (t18 * t46 + t19 * t44) * qJD(1)) + m(5) * (t101 * qJD(1) + t44 * t7 - t46 * t6); 0; 0; (-t68 * qJD(4) - t22 * t105 - t24 * t106) * t44 / 0.2e1 - (-t69 * qJD(4) + t23 * t105 + t25 * t106) * t46 / 0.2e1 + m(5) * ((-t44 * t6 - t46 * t7) * t35 - t101 * t31) + (t99 / 0.2e1 + t98 / 0.2e1) * t63 * qJD(4) + ((t88 / 0.2e1 + t90 / 0.2e1 - t17 * t97) * t46 + (t89 / 0.2e1 + t91 / 0.2e1 + t16 * t97) * t44) * qJD(1); m(5) * t1; 0; ((t26 * t44 + t27 * t46) * t1 + (t98 + t99) * t35 * t31) * t100 + (-t4 * t46 + t5 * t44) * t77 + t44 * ((t44 * t8 + (t4 + t103) * qJD(1)) * t44 + (t5 * qJD(1) + (t22 * t75 + t24 * t76) * t46 + (-t9 + (-t88 - t90) * qJD(4) + (t21 - t69) * qJD(1)) * t44) * t46) + (-t2 * t46 + t3 * t44) * t78 - t46 * ((t46 * t9 + (t3 + t102) * qJD(1)) * t46 + (t2 * qJD(1) + (-t23 * t75 - t25 * t76 + t79) * t44 + (-t8 + (t89 + t91) * qJD(4) - t68 * qJD(1)) * t46) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
