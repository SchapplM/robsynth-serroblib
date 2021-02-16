% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:48
% EndTime: 2021-01-15 19:14:54
% DurationCPUTime: 0.80s
% Computational Cost: add. (1079->131), mult. (2525->239), div. (0->0), fcn. (2331->6), ass. (0->76)
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t48 = sin(qJ(3));
t50 = cos(qJ(3));
t30 = t48 * t45 - t50 * t46;
t81 = pkin(6) + qJ(2);
t32 = t81 * t45;
t33 = t81 * t46;
t59 = -t50 * t32 - t48 * t33;
t53 = t30 * qJD(2) - t59 * qJD(3);
t31 = t50 * t45 + t48 * t46;
t37 = -t46 * pkin(2) - pkin(1);
t56 = -t30 * pkin(3) + t31 * pkin(7) - t37;
t95 = qJD(4) * t56 + t53;
t22 = -t48 * t32 + t50 * t33;
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t79 = t49 * t22 - t47 * t56;
t26 = t30 * qJD(3);
t27 = t31 * qJD(3);
t65 = t27 * pkin(3) + t26 * pkin(7);
t13 = t49 * t65;
t55 = -t26 * qJ(5) + qJD(4) * t22 + t31 * qJD(5);
t77 = qJ(5) * t31;
t71 = qJD(4) * t77;
t51 = -t55 * t49 + t13 + (t71 + t95) * t47;
t91 = t27 * pkin(4);
t1 = t51 + t91;
t74 = -t47 * t65 + t95 * t49;
t2 = t55 * t47 + t49 * t71 + t74;
t68 = -t47 * t22 - t49 * t56;
t5 = t30 * pkin(4) - t49 * t77 + t68;
t6 = -t47 * t77 + t79;
t94 = -t1 * t49 + t2 * t47 + (t47 * t5 - t49 * t6) * qJD(4);
t93 = 0.2e1 * t37;
t92 = 0.2e1 * qJD(4);
t90 = t49 * pkin(4);
t89 = t31 * t26;
t88 = t31 * t47;
t87 = t31 * t49;
t86 = t47 * t26;
t85 = t47 * t27;
t84 = t49 * t26;
t83 = t49 * t27;
t80 = -qJ(5) - pkin(7);
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t78 = t43 - t44;
t76 = qJD(4) * t47;
t40 = qJD(4) * t49;
t75 = -0.2e1 * pkin(3) * qJD(4);
t73 = pkin(4) * t76;
t72 = t47 * t40;
t38 = -pkin(3) - t90;
t70 = -t38 + t90;
t69 = -0.4e1 * t47 * t87;
t67 = t78 * qJD(4);
t66 = 0.2e1 * (t45 ^ 2 + t46 ^ 2) * qJD(2);
t64 = pkin(3) * t26 - pkin(7) * t27;
t63 = pkin(3) * t31 + pkin(7) * t30;
t60 = pkin(4) * t43 + t38 * t49;
t19 = t31 * t40 - t86;
t58 = t31 * t76 + t84;
t57 = t30 * t40 + t85;
t9 = t31 * qJD(2) + t22 * qJD(3);
t35 = t80 * t49;
t34 = t80 * t47;
t28 = t31 ^ 2;
t25 = -t47 * qJD(5) + t80 * t40;
t24 = -t49 * qJD(5) - t80 * t76;
t16 = -t30 * t76 + t83;
t10 = pkin(4) * t88 - t59;
t7 = t19 * pkin(4) + t9;
t4 = -t79 * qJD(4) + t47 * t53 + t13;
t3 = t22 * t76 + t74;
t8 = [0, 0, 0, 0, t66, qJ(2) * t66, -0.2e1 * t89, 0.2e1 * t26 * t30 - 0.2e1 * t31 * t27, 0, 0, 0, t27 * t93, -t26 * t93, -0.2e1 * t28 * t72 - 0.2e1 * t44 * t89, t78 * t28 * t92 - t26 * t69, -0.2e1 * t58 * t30 + 0.2e1 * t31 * t83, -0.2e1 * t19 * t30 - 0.2e1 * t31 * t85, 0.2e1 * t30 * t27, -0.2e1 * t19 * t59 + 0.2e1 * t68 * t27 + 0.2e1 * t4 * t30 + 0.2e1 * t9 * t88, -0.2e1 * t79 * t27 + 0.2e1 * t3 * t30 + 0.2e1 * t58 * t59 + 0.2e1 * t9 * t87, 0.2e1 * t1 * t30 + 0.2e1 * t19 * t10 + 0.2e1 * t5 * t27 + 0.2e1 * t7 * t88, -0.2e1 * t58 * t10 + 0.2e1 * t2 * t30 - 0.2e1 * t6 * t27 + 0.2e1 * t7 * t87, -0.2e1 * (-t47 * t6 - t49 * t5) * t26 + 0.2e1 * t94 * t31, 0.2e1 * t5 * t1 + 0.2e1 * t10 * t7 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, 0, 0, 0, 0, t16, -t57, t16, -t57, -(-t43 - t44) * t26, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, -t9, t53, -t31 * t67 - t47 * t84, qJD(4) * t69 + t78 * t26, t57, t16, 0, -t9 * t49 + t64 * t47 + (-t47 * t59 - t63 * t49) * qJD(4), t9 * t47 + t64 * t49 + (t63 * t47 - t49 * t59) * qJD(4), -t38 * t86 + t25 * t30 + t34 * t27 - t7 * t49 + (t10 * t47 + t60 * t31) * qJD(4), -t38 * t84 + t24 * t30 + t35 * t27 + t7 * t47 + (t10 * t49 + t70 * t88) * qJD(4), (-t25 * t31 + t26 * t34 - t2 + (t31 * t35 - t5) * qJD(4)) * t49 + (t24 * t31 - t26 * t35 - t1 + (t31 * t34 - t6) * qJD(4)) * t47, t1 * t34 + t10 * t73 + t2 * t35 - t6 * t24 + t5 * t25 + t7 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t24 + t49 * t25 + (-t34 * t47 - t35 * t49) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72, -0.2e1 * t67, 0, 0, 0, t47 * t75, t49 * t75, -0.2e1 * t70 * t76, t60 * t92, -0.2e1 * t24 * t49 - 0.2e1 * t25 * t47 + 0.2e1 * (-t34 * t49 + t35 * t47) * qJD(4), 0.2e1 * t35 * t24 + 0.2e1 * t34 * t25 + 0.2e1 * t38 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t19, t27, t4, t3, t51 + 0.2e1 * t91, t2, t58 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t40, -t76, -t40, 0, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t76, 0, -pkin(7) * t40, pkin(7) * t76, t25, t24, -pkin(4) * t40, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t58, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t40, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
