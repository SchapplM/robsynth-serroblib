% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:15
% DurationCPUTime: 0.44s
% Computational Cost: add. (596->89), mult. (1367->139), div. (0->0), fcn. (1163->6), ass. (0->66)
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t82 = t60 ^ 2 + t61 ^ 2;
t100 = t82 * qJD(3);
t101 = 0.2e1 * t100;
t64 = cos(qJ(2));
t81 = pkin(1) * qJD(2);
t79 = t64 * t81;
t48 = qJD(3) + t79;
t99 = t82 * t48;
t98 = 2 * qJD(5);
t97 = t64 * pkin(1);
t63 = sin(qJ(2));
t52 = t63 * pkin(1) + qJ(3);
t96 = -pkin(7) - t52;
t91 = cos(qJ(4));
t72 = qJD(4) * t91;
t62 = sin(qJ(4));
t80 = qJD(4) * t62;
t37 = t60 * t80 - t61 * t72;
t77 = t91 * t60;
t88 = t62 * t61;
t43 = t77 + t88;
t38 = t43 * qJD(4);
t19 = t38 * pkin(4) + t37 * qJ(5) - t43 * qJD(5);
t54 = t63 * t81;
t12 = t19 + t54;
t76 = t91 * t61;
t42 = t62 * t60 - t76;
t53 = -t61 * pkin(3) - pkin(2);
t28 = t42 * pkin(4) - t43 * qJ(5) + t53;
t25 = t28 - t97;
t95 = t12 * t42 + t25 * t38;
t94 = -t12 * t43 + t25 * t37;
t93 = t19 * t42 + t28 * t38;
t92 = -t19 * t43 + t28 * t37;
t90 = t53 * t37;
t89 = t53 * t38;
t87 = -pkin(7) - qJ(3);
t46 = t53 - t97;
t86 = t46 * t38 + t42 * t54;
t85 = -t46 * t37 + t43 * t54;
t78 = t62 * t96;
t75 = t87 * t60;
t57 = t61 * pkin(7);
t39 = t61 * t52 + t57;
t66 = t96 * t77;
t26 = t62 * t39 - t66;
t27 = t91 * t39 + t60 * t78;
t5 = -qJD(4) * t66 - t48 * t76 + (qJD(4) * t39 + t48 * t60) * t62;
t6 = t39 * t72 + t48 * t88 + (qJD(4) * t78 + t91 * t48) * t60;
t73 = -t26 * t37 - t27 * t38 + t5 * t42 + t6 * t43;
t71 = t91 * qJD(3);
t47 = t61 * qJ(3) + t57;
t65 = t91 * t75;
t21 = -qJD(4) * t65 - t61 * t71 + (qJD(3) * t60 + qJD(4) * t47) * t62;
t22 = t47 * t72 + qJD(3) * t88 + (t87 * t80 + t71) * t60;
t30 = t62 * t47 - t65;
t31 = t91 * t47 + t62 * t75;
t69 = t21 * t42 + t22 * t43 - t30 * t37 - t31 * t38;
t68 = t60 * t54;
t67 = t61 * t54;
t29 = -0.2e1 * t43 * t37;
t20 = pkin(4) * t37 - t38 * qJ(5) - t42 * qJD(5);
t9 = 0.2e1 * t37 * t42 - 0.2e1 * t43 * t38;
t1 = [0, 0, 0, 0, -0.2e1 * t54, -0.2e1 * t79, -0.2e1 * t67, 0.2e1 * t68, 0.2e1 * t99, 0.2e1 * (-pkin(2) - t97) * t54 + 0.2e1 * t52 * t99, t29, t9, 0, 0, 0, 0.2e1 * t86, 0.2e1 * t85, 0.2e1 * t95, 0.2e1 * t73, 0.2e1 * t94, 0.2e1 * t25 * t12 + 0.2e1 * t26 * t6 - 0.2e1 * t27 * t5; 0, 0, 0, 0, -t54, -t79, -t67, t68, t100 + t99, -pkin(2) * t54 + qJ(3) * t99 + t100 * t52, t29, t9, 0, 0, 0, t86 + t89, t85 - t90, t93 + t95, t69 + t73, t92 + t94, t12 * t28 + t25 * t19 - t27 * t21 + t26 * t22 + t6 * t30 - t5 * t31; 0, 0, 0, 0, 0, 0, 0, 0, t101, qJ(3) * t101, t29, t9, 0, 0, 0, 0.2e1 * t89, -0.2e1 * t90, 0.2e1 * t93, 0.2e1 * t69, 0.2e1 * t92, 0.2e1 * t28 * t19 - 0.2e1 * t31 * t21 + 0.2e1 * t30 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, t38, -t37, t38, 0, t37, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, t38, 0, t37, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, -t6, t5, -t6, t20, -t5, -t6 * pkin(4) - t5 * qJ(5) + t27 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, -t22, t21, -t22, t20, -t21, -t22 * pkin(4) - t21 * qJ(5) + t31 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, qJ(5) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
