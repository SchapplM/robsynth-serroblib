% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:06
% DurationCPUTime: 0.70s
% Computational Cost: add. (607->124), mult. (1328->228), div. (0->0), fcn. (917->4), ass. (0->85)
t95 = pkin(3) + pkin(6);
t42 = sin(qJ(2));
t85 = t42 * qJ(3);
t44 = cos(qJ(2));
t45 = -pkin(2) - pkin(7);
t90 = t44 * t45;
t99 = t85 - t90;
t21 = -pkin(1) - t99;
t29 = t95 * t42;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t98 = t43 * t21 + t41 * t29;
t40 = t44 ^ 2;
t58 = qJD(2) * (t42 ^ 2 - t40);
t37 = t41 ^ 2;
t88 = -t43 ^ 2 + t37;
t57 = t88 * qJD(4);
t77 = t42 * qJD(2);
t56 = pkin(2) * t77 - t42 * qJD(3);
t86 = qJ(3) * t44;
t11 = (pkin(7) * t42 - t86) * qJD(2) + t56;
t34 = t44 * qJD(2);
t32 = pkin(6) * t34;
t25 = pkin(3) * t34 + t32;
t5 = -qJD(4) * t98 - t41 * t11 + t25 * t43;
t97 = 0.2e1 * qJD(3);
t96 = 2 * qJD(5);
t51 = -pkin(4) * t41 + qJ(5) * t43;
t14 = qJD(4) * t51 + t41 * qJD(5);
t52 = pkin(4) * t43 + qJ(5) * t41;
t6 = t14 * t44 + (-t52 - t95) * t77;
t94 = t6 * t41;
t12 = qJD(4) * t52 - t43 * qJD(5) + qJD(3);
t93 = t12 * t44;
t24 = t95 * t77;
t92 = t24 * t41;
t91 = t42 * t45;
t30 = t95 * t44;
t84 = qJD(2) * t41;
t83 = qJD(2) * t43;
t82 = qJD(4) * t41;
t81 = qJD(4) * t43;
t80 = qJD(4) * t44;
t79 = qJD(4) * t45;
t78 = t30 * qJD(4);
t76 = t42 * qJD(5);
t75 = t44 * qJD(3);
t74 = qJ(3) * qJD(4);
t73 = qJ(5) * qJD(2);
t72 = qJD(2) * qJ(3);
t71 = -0.2e1 * pkin(1) * qJD(2);
t70 = pkin(4) * t34;
t69 = pkin(6) * t77;
t68 = t41 * t80;
t67 = t41 * t79;
t66 = t43 * t80;
t65 = t43 * t79;
t64 = t42 * t34;
t63 = t43 * t77;
t62 = t43 * t34;
t61 = t41 * t81;
t60 = t44 * t73;
t55 = t41 * t63;
t7 = qJ(5) * t42 + t98;
t49 = -t21 * t41 + t29 * t43;
t8 = -pkin(4) * t42 - t49;
t54 = t41 * t8 + t43 * t7;
t53 = -pkin(2) * t44 - t85;
t26 = qJ(3) - t51;
t48 = t26 * t44 + t91;
t47 = -t86 - t91;
t4 = -t11 * t43 + t21 * t82 - t25 * t41 - t29 * t81;
t46 = qJD(2) * t53 + t75;
t2 = -t4 + t60 + t76;
t3 = -t5 - t70;
t1 = qJD(4) * t54 + t2 * t41 - t3 * t43;
t31 = 0.2e1 * t64;
t28 = t45 * t62;
t27 = -pkin(1) + t53;
t18 = t41 * t77 - t66;
t17 = t34 * t41 + t42 * t81;
t16 = -t42 * t82 + t62;
t15 = -t44 * t72 + t56;
t10 = t44 * t52 + t30;
t9 = [0, 0, 0, t31, -0.2e1 * t58, 0, 0, 0, t42 * t71, t44 * t71, 0, 0.2e1 * t15 * t44 - 0.2e1 * t27 * t77, -0.2e1 * t15 * t42 - 0.2e1 * t27 * t34, 0.2e1 * t27 * t15, -0.2e1 * t37 * t64 + 0.2e1 * t40 * t61, -0.2e1 * t40 * t57 - 0.4e1 * t44 * t55, 0.2e1 * t41 * t58 - 0.2e1 * t42 * t66, 0.2e1 * t42 * t68 + 0.2e1 * t43 * t58, t31, 0.2e1 * (-t30 * t83 + t5) * t42 + 0.2e1 * (qJD(2) * t49 - t24 * t43 - t41 * t78) * t44, 0.2e1 * (t30 * t84 + t4) * t42 + 0.2e1 * (-qJD(2) * t98 - t43 * t78 + t92) * t44, 0.2e1 * (-t10 * t83 - t3) * t42 + 0.2e1 * (-qJD(2) * t8 - t10 * t82 + t6 * t43) * t44, 0.2e1 * t54 * t77 + 0.2e1 * (-t2 * t43 - t3 * t41 + (t41 * t7 - t43 * t8) * qJD(4)) * t44, 0.2e1 * (-t10 * t84 + t2) * t42 + 0.2e1 * (qJD(2) * t7 + t10 * t81 + t94) * t44, 0.2e1 * t10 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, t34, -t77, 0, -t32, t69, t46, t32, -t69, t46 * pkin(6), t44 * t57 + t55, 0.4e1 * t44 * t61 - t77 * t88, t16, -t17, 0, -t92 + t28 + (-t42 * t72 + t75) * t43 + (t30 * t43 + t41 * t47) * qJD(4), (qJD(4) * t47 - t24) * t43 + (qJD(2) * t99 - t75 - t78) * t41, t94 + t28 + (-t26 * t77 + t93) * t43 + (t10 * t43 - t41 * t48) * qJD(4), -t1, (qJD(4) * t48 - t6) * t43 + (qJD(4) * t10 + t93 + (-t26 * t42 + t90) * qJD(2)) * t41, t1 * t45 + t10 * t12 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(3) * t97, -0.2e1 * t61, 0.2e1 * t57, 0, 0, 0, 0.2e1 * qJD(3) * t41 + 0.2e1 * t43 * t74, 0.2e1 * qJD(3) * t43 - 0.2e1 * t41 * t74, 0.2e1 * t12 * t41 + 0.2e1 * t26 * t81, 0, -0.2e1 * t12 * t43 + 0.2e1 * t26 * t82, 0.2e1 * t26 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, t32, 0, 0, 0, 0, 0, t16, -t17, t16, 0, t17, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t63 + t68, t34, t5, t4, t5 + 0.2e1 * t70, (-pkin(4) * t77 + qJ(5) * t80) * t41 + (t42 * t73 + (pkin(4) * qJD(4) - qJD(5)) * t44) * t43, -t4 + 0.2e1 * t60 + 0.2e1 * t76, -pkin(4) * t3 + qJ(5) * t2 + qJD(5) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t81, 0, -t67, -t65, -t67, -t14, t65, t14 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t81, -t82, 0, t81, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, qJ(5) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t18, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
