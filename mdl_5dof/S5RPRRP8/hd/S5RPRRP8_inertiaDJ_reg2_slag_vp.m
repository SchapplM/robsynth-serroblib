% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:30
% EndTime: 2019-12-31 18:47:32
% DurationCPUTime: 0.75s
% Computational Cost: add. (409->88), mult. (858->126), div. (0->0), fcn. (535->4), ass. (0->59)
t77 = cos(qJ(3));
t81 = -pkin(1) - pkin(2);
t56 = t77 * t81;
t48 = sin(qJ(3));
t68 = qJD(3) * t48;
t16 = qJ(2) * t68 - t77 * qJD(2) - qJD(3) * t56;
t47 = sin(qJ(4));
t45 = t47 ^ 2;
t49 = cos(qJ(4));
t46 = t49 ^ 2;
t86 = t45 + t46;
t2 = t86 * t16;
t60 = qJD(3) * t77;
t88 = t86 * t60;
t85 = -t47 * pkin(4) + t49 * qJ(5);
t19 = -qJD(4) * t85 - t47 * qJD(5);
t53 = t49 * pkin(4) + t47 * qJ(5);
t18 = t53 * qJD(4) - t49 * qJD(5);
t84 = 0.2e1 * (-t45 + t46) * qJD(4);
t83 = 0.2e1 * qJD(2);
t82 = 0.2e1 * qJD(5);
t30 = t77 * qJ(2) + t48 * t81;
t28 = -pkin(7) + t30;
t80 = t2 * t28;
t17 = t48 * qJD(2) + t30 * qJD(3);
t5 = t17 - t19;
t78 = t19 - t5;
t76 = t17 * t47;
t75 = t17 * t49;
t72 = t2 * pkin(7);
t29 = -t48 * qJ(2) + t56;
t27 = pkin(3) - t29;
t13 = t27 + t53;
t31 = -pkin(3) - t53;
t71 = t13 - t31;
t70 = t88 * pkin(7);
t40 = t47 * qJD(4);
t41 = t49 * qJD(4);
t65 = -0.2e1 * pkin(3) * qJD(4);
t64 = pkin(7) * t40;
t63 = pkin(7) * t41;
t62 = t47 * t41;
t61 = t17 * t77;
t59 = qJD(4) * (pkin(3) + t27);
t58 = qJD(4) * t77;
t57 = -t2 * t48 + t88 * t28;
t51 = t86 * t77;
t36 = -0.2e1 * t62;
t35 = 0.2e1 * t62;
t24 = t47 * t58 + t49 * t68;
t23 = t47 * t68 - t49 * t58;
t22 = t48 * t41 + t47 * t60;
t21 = t51 * qJD(3);
t20 = t48 * t40 - t49 * t60;
t8 = 0.2e1 * (-t77 + t51) * t68;
t4 = -t47 * t16 + t28 * t41;
t3 = t49 * t16 + t28 * t40;
t1 = -0.2e1 * t2;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, qJ(2) * t83, 0, 0, 0, 0, 0, 0, 0.2e1 * t17, -0.2e1 * t16, 0, -0.2e1 * t30 * t16 - 0.2e1 * t29 * t17, t35, t84, 0, t36, 0, 0, -0.2e1 * t27 * t40 + 0.2e1 * t75, -0.2e1 * t27 * t41 - 0.2e1 * t76, -t1, 0.2e1 * t27 * t17 - 0.2e1 * t80, t35, 0, -t84, 0, 0, t36, -0.2e1 * t13 * t40 + 0.2e1 * t5 * t49, -t1, 0.2e1 * t13 * t41 + 0.2e1 * t5 * t47, 0.2e1 * t13 * t5 - 0.2e1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t60, 0, -t61 - t16 * t48 + (-t29 * t48 + t77 * t30) * qJD(3), 0, 0, 0, 0, 0, 0, t24, -t23, -t21, t27 * t68 + t57 - t61, 0, 0, 0, 0, 0, 0, t24, -t21, t23, t13 * t68 - t5 * t77 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t36, -t84, 0, t35, 0, 0, t47 * t59 - t75, t49 * t59 + t76, -t2, -t17 * pkin(3) - t72, t36, 0, t84, 0, 0, t35, t71 * t40 + t78 * t49, -t2, -t71 * t41 + t78 * t47, t13 * t19 + t5 * t31 - t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t60, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t23, t21, -pkin(3) * t68 + t70, 0, 0, 0, 0, 0, 0, -t24, t21, -t23, -t77 * t19 + t31 * t68 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t84, 0, t36, 0, 0, t47 * t65, t49 * t65, 0, 0, t35, 0, -t84, 0, 0, t36, -0.2e1 * t19 * t49 + 0.2e1 * t31 * t40, 0, -0.2e1 * t19 * t47 - 0.2e1 * t31 * t41, 0.2e1 * t31 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, t40, 0, -t4, t3, 0, 0, 0, -t41, 0, 0, -t40, 0, -t4, t18, -t3, -t16 * t85 - t18 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t20, -t18 * t48 + t85 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t40, 0, -t63, t64, 0, 0, 0, t41, 0, 0, t40, 0, -t63, -t18, -t64, -t18 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, qJ(5) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
