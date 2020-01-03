% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR6
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
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:43
% EndTime: 2019-12-31 16:52:44
% DurationCPUTime: 0.51s
% Computational Cost: add. (754->65), mult. (1551->99), div. (0->0), fcn. (1792->6), ass. (0->65)
t118 = qJD(3) + qJD(4);
t103 = cos(qJ(3));
t53 = sin(pkin(7));
t54 = cos(pkin(7));
t56 = sin(qJ(3));
t42 = -t103 * t54 + t56 * t53;
t95 = pkin(5) + qJ(2);
t46 = t95 * t53;
t47 = t95 * t54;
t58 = -t103 * t47 + t56 * t46;
t24 = -t42 * pkin(6) - t58;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t44 = t103 * t53 + t56 * t54;
t63 = t103 * t46 + t56 * t47;
t59 = -t44 * pkin(6) - t63;
t117 = t118 * (-t57 * t24 - t55 * t59);
t116 = t118 * (t55 * t24 - t57 * t59);
t39 = t57 * t44;
t99 = t55 * t42;
t105 = t39 - t99;
t30 = t57 * t42 + t55 * t44;
t108 = -t105 ^ 2 + t30 ^ 2;
t115 = t108 * qJD(1);
t79 = t30 * qJD(4);
t9 = -t30 * qJD(3) - t79;
t110 = qJD(1) * t30;
t109 = qJD(2) * t30;
t80 = t105 * qJD(1);
t65 = t39 / 0.2e1;
t104 = pkin(3) * t44;
t48 = t53 ^ 2 + t54 ^ 2;
t94 = pkin(3) * qJD(4);
t93 = qJD(3) * pkin(3);
t50 = -t54 * pkin(2) - pkin(1);
t34 = t42 * pkin(3) + t50;
t6 = -t30 * t104 - t105 * t34;
t89 = t6 * qJD(1);
t7 = -t104 * t105 + t30 * t34;
t88 = t7 * qJD(1);
t86 = qJD(1) * t34;
t13 = 0.2e1 * t65 - t99;
t84 = t13 * qJD(1);
t25 = t42 ^ 2 - t44 ^ 2;
t83 = t25 * qJD(1);
t28 = t65 - t39 / 0.2e1;
t82 = t28 * qJD(1);
t81 = t28 * qJD(4);
t76 = t105 * qJD(4);
t75 = t42 * qJD(1);
t40 = t42 * qJD(3);
t74 = t44 * qJD(1);
t73 = t44 * qJD(3);
t45 = t48 * qJ(2);
t72 = t45 * qJD(1);
t71 = t48 * qJD(1);
t70 = t30 * t80;
t69 = t105 * t110;
t68 = t30 * t86;
t67 = t105 * t86;
t66 = t42 * t74;
t62 = pkin(3) * t118;
t61 = qJD(1) * t50 + qJD(2);
t60 = qJD(3) * t105 + t13 * qJD(4);
t1 = [0, 0, 0, 0, 0, t48 * qJD(2), t45 * qJD(2), -t42 * t73, t25 * qJD(3), 0, 0, 0, t50 * t73, -t50 * t40, t9 * t105, t118 * t108, 0, 0, 0, -t6 * qJD(3) + t34 * t76, -t7 * qJD(3) - t34 * t79; 0, 0, 0, 0, 0, t71, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0; 0, 0, 0, 0, 0, 0, 0, -t66, t83, -t40, -t73, 0, t58 * qJD(3) + t50 * t74, t63 * qJD(3) - t50 * t75, -t69, t115, t9, -t60, 0, -t89 + t117, -t88 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t115, t9, -t13 * qJD(3) - t76, 0, t28 * qJD(2) + t117 + t67, -t68 + t116; 0, 0, 0, 0, 0, -t71, -t72, 0, 0, 0, 0, 0, t73, -t40, 0, 0, 0, 0, 0, t60, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t75, 0, 0, 0, 0, 0, t80, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t110; 0, 0, 0, 0, 0, 0, 0, t66, -t83, 0, 0, 0, -t61 * t44, t61 * t42, t69, -t115, 0, -t81, 0, -qJD(2) * t105 + t89, t88 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t75, 0, 0, 0, 0, 0, -t80, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t94, -t57 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, -t55 * t62, -t57 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t115, 0, t28 * qJD(3), 0, -t13 * qJD(2) - t67, t68 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, t55 * t93, t57 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
