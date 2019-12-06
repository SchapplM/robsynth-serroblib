% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:15
% EndTime: 2019-12-05 15:45:17
% DurationCPUTime: 0.45s
% Computational Cost: add. (573->69), mult. (1250->109), div. (0->0), fcn. (1404->8), ass. (0->66)
t66 = sin(pkin(9));
t67 = cos(pkin(9));
t91 = sin(qJ(2));
t92 = cos(qJ(2));
t54 = -t66 * t91 + t67 * t92;
t55 = -t66 * t92 - t67 * t91;
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t78 = t69 * t54 - t71 * t55;
t84 = t78 * qJD(4);
t85 = t78 * qJD(2);
t108 = -t85 - t84;
t81 = qJD(2) + qJD(4);
t52 = t71 * t54;
t90 = t69 * t55;
t33 = t52 + t90;
t107 = 0.2e1 * t33;
t70 = cos(qJ(5));
t95 = -t70 / 0.2e1;
t102 = t95 * t107;
t106 = qJD(5) * t102;
t68 = sin(qJ(5));
t97 = -t68 / 0.2e1;
t103 = t97 * t107;
t105 = qJD(5) * t103;
t104 = t81 * t33;
t101 = 0.2e1 * t78;
t100 = t70 * t81;
t58 = -t68 ^ 2 + t70 ^ 2;
t99 = t81 * t58;
t62 = t67 * pkin(2) + pkin(3);
t93 = pkin(2) * t66;
t50 = -t71 * t62 + t69 * t93;
t44 = -pkin(4) + t50;
t98 = t44 + t50;
t96 = t68 / 0.2e1;
t89 = pkin(4) * qJD(4);
t87 = qJD(2) * t44;
t86 = qJD(5) * t68;
t65 = qJD(5) * t70;
t83 = t50 * qJD(2);
t51 = t69 * t62 + t71 * t93;
t82 = t51 * qJD(2);
t43 = t51 * qJD(4);
t79 = t50 / 0.2e1 + pkin(4) / 0.2e1 - t44 / 0.2e1;
t77 = -t33 / 0.2e1 + t90 / 0.2e1 + t52 / 0.2e1;
t16 = t79 * t68;
t76 = t16 * qJD(2) + t68 * t89;
t17 = t79 * t70;
t75 = t17 * qJD(2) + t70 * t89;
t2 = t77 * t68;
t74 = t2 * qJD(1) - t68 * t87;
t73 = t68 * t82;
t6 = t77 * t70;
t72 = t6 * qJD(1) - t70 * t87;
t59 = t68 * t65;
t57 = t58 * qJD(5);
t53 = t68 * t100;
t45 = pkin(7) + t51;
t42 = t50 * qJD(4);
t39 = t68 * t43;
t19 = pkin(4) * t95 + t98 * t70 / 0.2e1;
t18 = pkin(4) * t97 + t98 * t96;
t8 = t95 * t101;
t3 = t96 * t101;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t91, -qJD(2) * t92, (t54 * t66 + t55 * t67) * qJD(2) * pkin(2), 0, t108, -t104, 0, 0, 0, 0, 0, t8 * qJD(4) - t70 * t85 + t105, t3 * qJD(4) + t68 * t85 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t108, -t104, 0, 0, 0, 0, 0, t8 * qJD(2) - t70 * t84 + t105, t3 * qJD(2) + t68 * t84 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t103 - t65 * t78, t81 * t102 + t78 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(5), -t6 * qJD(5); 0, 0, 0, 0, 0, 0, -t43, t42, t59, t57, 0, 0, 0, -t70 * t43 + t44 * t86, t44 * t65 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t43 - t82, t42 + t83, t59, t57, 0, 0, 0, t18 * qJD(5) - t51 * t100, t19 * qJD(5) + t39 + t73; 0, 0, 0, 0, 0, 0, 0, 0, t53, t99, t65, -t86, 0, t18 * qJD(4) - t45 * t65 - t74, t19 * qJD(4) + t45 * t86 - t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t82, -t83, t59, t57, 0, 0, 0, -t16 * qJD(5) + t70 * t82, -t17 * qJD(5) - t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t59, t57, 0, 0, 0, -pkin(4) * t86, -pkin(4) * t65; 0, 0, 0, 0, 0, 0, 0, 0, t53, t99, t65, -t86, 0, -pkin(7) * t65 - t76, pkin(7) * t86 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(2), t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t99, 0, 0, 0, t16 * qJD(4) + t74, t17 * qJD(4) + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t99, 0, 0, 0, t76, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
