% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP10
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
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:45
% EndTime: 2021-01-15 21:02:50
% DurationCPUTime: 0.80s
% Computational Cost: add. (635->134), mult. (1409->247), div. (0->0), fcn. (980->4), ass. (0->87)
t99 = pkin(3) + pkin(6);
t100 = pkin(2) + pkin(7);
t89 = qJ(5) + t100;
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t92 = qJ(3) * t47;
t104 = t100 * t45 - t92;
t90 = t45 * qJ(3);
t103 = t100 * t47 + t90;
t43 = t47 ^ 2;
t63 = qJD(2) * (t45 ^ 2 - t43);
t44 = sin(qJ(4));
t40 = t44 ^ 2;
t46 = cos(qJ(4));
t95 = -t46 ^ 2 + t40;
t62 = t95 * qJD(4);
t77 = t47 * qJD(3);
t29 = t99 * t47;
t82 = t29 * qJD(4);
t102 = t103 * qJD(2) - t77 - t82;
t101 = 0.2e1 * qJD(3);
t98 = pkin(4) * t46;
t97 = t44 * t45;
t21 = -pkin(1) - t103;
t71 = t99 * t45;
t96 = t46 * t21 + t44 * t71;
t93 = qJ(3) * t44;
t91 = qJ(5) * t47;
t88 = qJD(2) * t44;
t87 = qJD(2) * t46;
t86 = qJD(4) * t44;
t85 = qJD(4) * t46;
t84 = qJD(4) * t47;
t83 = qJD(4) * t100;
t81 = t44 * qJD(5);
t80 = t45 * qJD(2);
t79 = t45 * qJD(3);
t78 = t46 * qJD(5);
t37 = t47 * qJD(2);
t76 = qJ(3) * qJD(4);
t75 = -0.2e1 * pkin(1) * qJD(2);
t74 = pkin(4) * t86;
t73 = pkin(6) * t80;
t72 = pkin(6) * t37;
t70 = t44 * t84;
t69 = t46 * t84;
t68 = t44 * t37;
t67 = t45 * t37;
t66 = t46 * t80;
t65 = t44 * t85;
t64 = qJ(5) * t84;
t27 = t89 * t46;
t24 = t46 * t71;
t61 = t44 * t66;
t60 = t99 * t37;
t36 = t44 * pkin(4) + qJ(3);
t8 = -pkin(4) * t70 + (-t98 - t99) * t80;
t59 = -t36 * t84 + t8;
t58 = -t21 * t85 + t46 * t60;
t6 = t45 * pkin(4) + t24 + (-t21 + t91) * t44;
t7 = -t46 * t91 + t96;
t57 = -t44 * t6 + t46 * t7;
t56 = -t47 * pkin(2) - t90;
t54 = -qJD(4) * t99 + qJD(3);
t4 = t21 * t86 - t44 * t60 - qJD(4) * t24 - t46 * (t104 * qJD(2) - t79);
t18 = t44 * t80 - t69;
t25 = t99 * t80;
t52 = t104 * qJD(4) - t25;
t14 = t47 * t98 + t29;
t34 = pkin(4) * t85 + qJD(3);
t51 = -qJD(4) * t14 - t34 * t47 + t36 * t80;
t50 = t56 * qJD(2) + t77;
t49 = t46 * t64 + t47 * t81 + (-t89 * qJD(2) + t54) * t97 + t58;
t33 = 0.2e1 * t67;
t28 = -pkin(1) + t56;
t26 = t89 * t44;
t17 = -t45 * t85 - t68;
t16 = t66 + t70;
t15 = t46 * t37 - t45 * t86;
t13 = -t79 + (pkin(2) * t45 - t92) * qJD(2);
t11 = -qJD(4) * t27 - t81;
t10 = t89 * t86 - t78;
t5 = qJ(3) * t68 + (-t100 * qJD(2) + t54) * t97 + t58;
t3 = t10 * t46 + t11 * t44 + (-t26 * t46 + t27 * t44) * qJD(4);
t2 = -qJ(5) * t66 - t44 * t64 + t47 * t78 + t4;
t1 = (pkin(4) + t93) * t37 + t49;
t9 = [0, 0, 0, t33, -0.2e1 * t63, 0, 0, 0, t45 * t75, t47 * t75, 0, 0.2e1 * t13 * t47 - 0.2e1 * t28 * t80, -0.2e1 * t13 * t45 - 0.2e1 * t28 * t37, 0.2e1 * t28 * t13, -0.2e1 * t40 * t67 + 0.2e1 * t43 * t65, -0.2e1 * t43 * t62 - 0.4e1 * t47 * t61, 0.2e1 * t44 * t63 - 0.2e1 * t45 * t69, 0.2e1 * t45 * t70 + 0.2e1 * t46 * t63, t33, 0.2e1 * (-t29 * t87 + t5) * t45 + 0.2e1 * ((-t44 * t21 + t24) * qJD(2) - t25 * t46 - t44 * t82) * t47, 0.2e1 * (t29 * t88 + t4) * t45 + 0.2e1 * (-t96 * qJD(2) + t25 * t44 - t46 * t82) * t47, 0.2e1 * (-t14 * t87 + t1) * t45 + 0.2e1 * (qJD(2) * t6 - t14 * t86 + t8 * t46) * t47, 0.2e1 * (t14 * t88 + t2) * t45 + 0.2e1 * (-qJD(2) * t7 - t14 * t85 - t8 * t44) * t47, 0.2e1 * t57 * t80 + 0.2e1 * (t1 * t44 + t2 * t46 + (t44 * t7 + t46 * t6) * qJD(4)) * t47, 0.2e1 * t6 * t1 + 0.2e1 * t14 * t8 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, t37, -t80, 0, -t72, t73, t50, t72, -t73, t50 * pkin(6), t47 * t62 + t61, 0.4e1 * t47 * t65 - t95 * t80, t15, t17, 0, -t102 * t46 + t52 * t44, t102 * t44 + t52 * t46, t10 * t45 - t27 * t37 + t59 * t44 - t51 * t46, -t11 * t45 + t26 * t37 + t51 * t44 + t59 * t46, (-t26 * t80 - t11 * t47 - t1 + (-t27 * t47 - t7) * qJD(4)) * t46 + (t27 * t80 + t10 * t47 + t2 + (-t26 * t47 + t6) * qJD(4)) * t44, -t1 * t27 + t6 * t10 + t7 * t11 + t14 * t34 + t2 * t26 + t8 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, qJ(3) * t101, -0.2e1 * t65, 0.2e1 * t62, 0, 0, 0, 0.2e1 * qJD(3) * t44 + 0.2e1 * t46 * t76, 0.2e1 * qJD(3) * t46 - 0.2e1 * t44 * t76, 0.2e1 * t34 * t44 + 0.2e1 * t36 * t85, 0.2e1 * t34 * t46 - 0.2e1 * t36 * t86, -0.2e1 * t3, -0.2e1 * t27 * t10 - 0.2e1 * t26 * t11 + 0.2e1 * t36 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, t72, 0, 0, 0, 0, 0, t15, t17, t15, t17, 0, t57 * qJD(4) + t1 * t46 - t2 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t16, t37, t5, t4, (0.2e1 * pkin(4) + t93) * t37 + t49, t2, -t18 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, t44 * t83, t46 * t83, t10, -t11, t74, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, -t86, -t85, 0, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t18, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t86, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
