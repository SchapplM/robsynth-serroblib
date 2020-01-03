% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:37
% EndTime: 2019-12-31 16:50:38
% DurationCPUTime: 0.47s
% Computational Cost: add. (264->86), mult. (696->150), div. (0->0), fcn. (543->6), ass. (0->84)
t43 = sin(qJ(4));
t38 = t43 ^ 2;
t45 = cos(qJ(4));
t40 = t45 ^ 2;
t29 = t40 - t38;
t44 = sin(qJ(3));
t75 = t44 * qJD(1);
t70 = t45 * t75;
t98 = t29 * qJD(3) - 0.2e1 * t43 * t70;
t46 = cos(qJ(3));
t95 = t46 * pkin(6);
t96 = t44 * pkin(3);
t24 = -t95 + t96;
t97 = -t24 / 0.2e1;
t34 = sin(pkin(7)) * pkin(1) + pkin(5);
t94 = t34 * t43;
t93 = t34 * t45;
t92 = t43 * t46;
t91 = t45 * t24;
t39 = t44 ^ 2;
t90 = t45 * t39;
t41 = t46 ^ 2;
t30 = t41 - t39;
t23 = t44 * t94;
t35 = -cos(pkin(7)) * pkin(1) - pkin(2);
t57 = -t46 * pkin(3) - t44 * pkin(6);
t14 = t35 + t57;
t72 = t34 * t92;
t8 = -t45 * t14 + t72;
t1 = t8 * t44 + (-t23 + t91) * t46;
t89 = t1 * qJD(1);
t71 = t46 * t93;
t9 = t43 * t14 + t71;
t2 = t24 * t92 + (-t9 + t71) * t44;
t88 = t2 * qJD(1);
t3 = -t39 * t94 - t8 * t46;
t87 = t3 * qJD(1);
t4 = -t34 * t90 - t9 * t46;
t86 = t4 * qJD(1);
t85 = qJD(3) * t43;
t84 = qJD(3) * t44;
t83 = qJD(3) * t45;
t82 = qJD(3) * t46;
t81 = qJD(4) * t43;
t80 = qJD(4) * t45;
t79 = qJD(4) * t46;
t20 = t30 * t43;
t78 = t20 * qJD(1);
t21 = t45 * t41 - t90;
t77 = t21 * qJD(1);
t76 = t30 * qJD(1);
t74 = t44 * qJD(4);
t73 = t46 * qJD(1);
t69 = t43 * t79;
t68 = t45 * t79;
t67 = t35 * t75;
t66 = t35 * t73;
t65 = t43 * t80;
t64 = t43 * t83;
t63 = t44 * t82;
t62 = t44 * t73;
t61 = t44 * t83;
t59 = t43 * t61;
t58 = -qJD(4) + t73;
t56 = t58 * t44;
t55 = t95 / 0.2e1 - t96 / 0.2e1;
t50 = t97 + t55;
t10 = t50 * t43;
t54 = pkin(3) * t83 + t10 * qJD(1);
t11 = t50 * t45;
t53 = pkin(3) * t85 - t11 * qJD(1);
t52 = t45 * t56;
t15 = (t38 / 0.2e1 - t40 / 0.2e1) * t44;
t51 = -t15 * qJD(1) + t64;
t49 = t43 * qJD(1) * t90 + t15 * qJD(3);
t19 = t29 * t39;
t48 = t19 * qJD(1) + 0.2e1 * t59;
t36 = t84 / 0.2e1;
t33 = t43 * t84;
t18 = (t73 - qJD(4) / 0.2e1) * t44;
t12 = t15 * qJD(4);
t6 = t23 + t91 / 0.2e1 + t55 * t45;
t5 = t44 * t93 + (-t55 + t97) * t43;
t7 = [0, 0, 0, 0, t63, t30 * qJD(3), 0, 0, 0, t35 * t84, t35 * t82, -t39 * t65 + t40 * t63, -t19 * qJD(4) - 0.2e1 * t46 * t59, -t21 * qJD(3) + t44 * t69, t20 * qJD(3) + t44 * t68, -t63, -t1 * qJD(3) - t4 * qJD(4), t2 * qJD(3) + t3 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t62, t76, t82, -t84, 0, -t34 * t82 + t67, t34 * t84 + t66, -t12 + (t40 * t75 + t64) * t46, -0.2e1 * t44 * t65 + t98 * t46, t33 - t77, t61 + t78, -t18, -t89 + (t43 * t57 - t71) * qJD(3) + t6 * qJD(4), t88 + (t45 * t57 + t72) * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, t43 * t56, t52, t36, t6 * qJD(3) - t9 * qJD(4) - t86, t5 * qJD(3) + t8 * qJD(4) + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t82, 0, 0, 0, 0, 0, -t61 - t69, t33 - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t82 - t45 * t74, t43 * t74 - t45 * t82; 0, 0, 0, 0, -t62, -t76, 0, 0, 0, -t67, -t66, -t40 * t62 - t12, 0.2e1 * t43 * t52, -t68 + t77, t69 - t78, t18, t11 * qJD(4) + t89, -t10 * qJD(4) - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t29 * qJD(4), 0, 0, 0, -pkin(3) * t81, -pkin(3) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t98, -t58 * t45, t58 * t43, -t75 / 0.2e1, -pkin(6) * t80 - t53, pkin(6) * t81 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, (-t43 * t75 + t83) * t46, (-t70 - t85) * t46, t36, -t11 * qJD(3) + t86, t10 * qJD(3) - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t98, t45 * t73, -t43 * t73, t75 / 0.2e1, t53, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
