% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP1
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 18:00:01
% DurationCPUTime: 0.66s
% Computational Cost: add. (773->90), mult. (1566->137), div. (0->0), fcn. (1334->4), ass. (0->59)
t45 = sin(qJ(3));
t46 = cos(qJ(4));
t47 = cos(qJ(3));
t75 = sin(qJ(4));
t31 = t46 * t45 + t75 * t47;
t79 = qJD(3) + qJD(4);
t16 = t79 * t31;
t41 = t46 * pkin(3) + pkin(4);
t61 = t75 * t45;
t72 = t46 * t47;
t32 = -t61 + t72;
t59 = qJD(4) * t75;
t54 = t32 * t59;
t17 = -qJD(3) * t61 - t45 * t59 + t72 * t79;
t63 = t75 * pkin(3);
t69 = qJD(4) * t46;
t64 = pkin(3) * t69;
t70 = t17 * t63 + t31 * t64;
t83 = -pkin(3) * t54 - t41 * t16 + t70;
t82 = t70 - (t16 * t46 + t54) * pkin(3);
t73 = t32 * t16;
t74 = t31 * t17;
t81 = 0.2e1 * t73 - 0.2e1 * t74;
t78 = 2 * qJD(2);
t77 = -pkin(1) - pkin(6);
t76 = t16 * pkin(4);
t65 = pkin(7) - t77;
t58 = t65 * t47;
t28 = t75 * t58;
t33 = t65 * t45;
t19 = -t46 * t33 - t28;
t39 = t45 * pkin(3) + qJ(2);
t68 = t45 * qJD(3);
t67 = t47 * qJD(3);
t34 = pkin(3) * t67 + qJD(2);
t66 = qJ(2) * qJD(3);
t62 = t45 * t67;
t60 = t77 * qJD(3);
t55 = pkin(3) * t59;
t30 = t46 * t58;
t18 = t75 * t33 - t30;
t52 = t65 * t68;
t6 = t79 * t30 - t33 * t59 - t75 * t52;
t10 = -t32 * qJ(5) + t18;
t11 = -t31 * qJ(5) + t19;
t2 = t17 * qJ(5) + t31 * qJD(5) + t6;
t7 = t28 * t79 + t33 * t69 + t46 * t52;
t3 = t16 * qJ(5) - t32 * qJD(5) + t7;
t49 = t10 * t16 - t11 * t17 + t2 * t31 - t3 * t32;
t48 = t18 * t16 - t19 * t17 + t6 * t31 - t7 * t32;
t44 = qJ(2) * t78;
t38 = -0.2e1 * t64;
t37 = -0.2e1 * t55;
t20 = t31 * pkin(4) + t39;
t12 = t17 * pkin(4) + t34;
t9 = -0.2e1 * t73;
t8 = 0.2e1 * t74;
t4 = 0.2e1 * t16 * t31 - 0.2e1 * t32 * t17;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t44, -0.2e1 * t62, 0.2e1 * (t45 ^ 2 - t47 ^ 2) * qJD(3), 0, 0.2e1 * t62, 0, 0, 0.2e1 * qJD(2) * t45 + 0.2e1 * t47 * t66, 0.2e1 * qJD(2) * t47 - 0.2e1 * t45 * t66, 0, t44, t9, t4, 0, t8, 0, 0, 0.2e1 * t39 * t17 + 0.2e1 * t34 * t31, -0.2e1 * t39 * t16 + 0.2e1 * t34 * t32, 0.2e1 * t48, 0.2e1 * t18 * t7 - 0.2e1 * t19 * t6 + 0.2e1 * t39 * t34, t9, t4, 0, t8, 0, 0, 0.2e1 * t12 * t31 + 0.2e1 * t20 * t17, 0.2e1 * t12 * t32 - 0.2e1 * t20 * t16, 0.2e1 * t49, 0.2e1 * t10 * t3 - 0.2e1 * t11 * t2 + 0.2e1 * t20 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t48, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, -t67, 0, -t45 * t60, -t47 * t60, 0, 0, 0, 0, -t16, 0, -t17, 0, t7, t6, -t82, (-t75 * t6 + t46 * t7 + (-t75 * t18 + t19 * t46) * qJD(4)) * pkin(3), 0, 0, -t16, 0, -t17, 0, t3, t2, -t83, t3 * t41 + (-t75 * t2 + (-t75 * t10 + t11 * t46) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, t82, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t38, 0, 0, 0, 0, 0, 0, 0, 0, t37, t38, 0, 0.2e1 * (-t75 * t41 + t46 * t63) * qJD(4) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, -t17, 0, t7, t6, 0, 0, 0, 0, -t16, 0, -t17, 0, t3, t2, t76, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t64, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t64, 0, -pkin(4) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
