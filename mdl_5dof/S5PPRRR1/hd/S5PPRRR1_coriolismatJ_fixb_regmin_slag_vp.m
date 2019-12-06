% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:53
% EndTime: 2019-12-05 15:12:54
% DurationCPUTime: 0.44s
% Computational Cost: add. (460->67), mult. (1068->114), div. (0->0), fcn. (1208->8), ass. (0->66)
t77 = qJD(3) + qJD(4);
t58 = sin(pkin(9));
t61 = sin(qJ(3));
t83 = cos(pkin(9));
t89 = cos(qJ(3));
t42 = t61 * t58 - t89 * t83;
t43 = t89 * t58 + t61 * t83;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t14 = t60 * t42 - t63 * t43;
t78 = t14 * qJD(4);
t79 = t14 * qJD(3);
t98 = t79 + t78;
t40 = t63 * t42;
t88 = t60 * t43;
t13 = t40 + t88;
t97 = t77 * t13;
t62 = cos(qJ(5));
t91 = t62 / 0.2e1;
t24 = t13 * t91;
t59 = sin(qJ(5));
t93 = t59 / 0.2e1;
t23 = t13 * t93;
t50 = -t59 ^ 2 + t62 ^ 2;
t96 = t77 * t50;
t95 = -0.2e1 * t14;
t94 = -t59 / 0.2e1;
t92 = -t62 / 0.2e1;
t90 = t63 * pkin(3);
t87 = pkin(3) * qJD(4);
t86 = pkin(4) * qJD(4);
t85 = qJD(3) * pkin(3);
t54 = -pkin(4) - t90;
t82 = qJD(3) * t54;
t81 = qJD(5) * t59;
t57 = qJD(5) * t62;
t76 = t60 * t87;
t75 = t60 * t85;
t74 = -t90 / 0.2e1;
t72 = pkin(3) * t77;
t71 = t60 * t72;
t70 = t74 + pkin(4) / 0.2e1 - t54 / 0.2e1;
t69 = t13 / 0.2e1 - t88 / 0.2e1 - t40 / 0.2e1;
t30 = t70 * t59;
t68 = t30 * qJD(3) + t59 * t86;
t31 = t70 * t62;
t67 = t31 * qJD(3) + t62 * t86;
t2 = t69 * t59;
t66 = t2 * qJD(1) - t59 * t82;
t6 = t69 * t62;
t65 = t6 * qJD(1) - t62 * t82;
t64 = t59 * t75;
t53 = t60 * pkin(3) + pkin(7);
t51 = t59 * t57;
t49 = t59 * t76;
t46 = t50 * qJD(5);
t41 = t77 * t62 * t59;
t33 = pkin(4) * t92 + t54 * t91 + t62 * t74;
t32 = pkin(4) * t94 + t54 * t93 + t59 * t74;
t12 = -t13 * t92 + t24;
t11 = 0.2e1 * t24;
t10 = -t13 * t94 + t23;
t9 = 0.2e1 * t23;
t8 = t95 * t92;
t3 = t95 * t93;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t43 * qJD(3), t42 * qJD(3), 0, t98, t97, 0, 0, 0, 0, 0, t8 * qJD(4) + t10 * qJD(5) + t62 * t79, t3 * qJD(4) + t12 * qJD(5) - t59 * t79; 0, 0, 0, 0, 0, 0, t98, t97, 0, 0, 0, 0, 0, t8 * qJD(3) + t9 * qJD(5) + t62 * t78, t3 * qJD(3) + t11 * qJD(5) - t59 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(3) + t9 * qJD(4) + t14 * t57, t12 * qJD(3) + t11 * qJD(4) - t14 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(5), -t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t76, -t63 * t87, t51, t46, 0, 0, 0, t54 * t81 - t62 * t76, t54 * t57 + t49; 0, 0, 0, 0, 0, 0, -t71, -t63 * t72, t51, t46, 0, 0, 0, t32 * qJD(5) - t62 * t71, t33 * qJD(5) + t49 + t64; 0, 0, 0, 0, 0, 0, 0, 0, t41, t96, t57, -t81, 0, t32 * qJD(4) - t53 * t57 - t66, t33 * qJD(4) + t53 * t81 - t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t75, t63 * t85, t51, t46, 0, 0, 0, -t30 * qJD(5) + t62 * t75, -t31 * qJD(5) - t64; 0, 0, 0, 0, 0, 0, 0, 0, t51, t46, 0, 0, 0, -pkin(4) * t81, -pkin(4) * t57; 0, 0, 0, 0, 0, 0, 0, 0, t41, t96, t57, -t81, 0, -pkin(7) * t57 - t68, pkin(7) * t81 - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(3), t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t96, 0, 0, 0, t30 * qJD(4) + t66, t31 * qJD(4) + t65; 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t96, 0, 0, 0, t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
