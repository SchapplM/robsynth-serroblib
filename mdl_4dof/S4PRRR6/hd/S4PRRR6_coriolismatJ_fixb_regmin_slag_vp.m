% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:07
% EndTime: 2019-12-31 16:35:08
% DurationCPUTime: 0.34s
% Computational Cost: add. (285->61), mult. (725->99), div. (0->0), fcn. (720->6), ass. (0->62)
t65 = qJD(3) + qJD(4);
t48 = sin(qJ(3));
t89 = pkin(5) + pkin(6);
t40 = t89 * t48;
t50 = cos(qJ(3));
t41 = t89 * t50;
t47 = sin(qJ(4));
t51 = cos(qJ(2));
t86 = cos(qJ(4));
t59 = t86 * t48;
t84 = t47 * t50;
t53 = -t84 / 0.2e1 - t59 / 0.2e1;
t54 = t59 + t84;
t16 = (t54 / 0.2e1 + t53) * t51;
t68 = t16 * qJD(1);
t92 = -t68 + t65 * (t47 * t40 - t86 * t41);
t58 = t86 * t50;
t85 = t47 * t48;
t33 = -t58 + t85;
t52 = -t58 / 0.2e1 + t85 / 0.2e1;
t17 = (-t33 / 0.2e1 + t52) * t51;
t67 = t17 * qJD(1);
t91 = -t67 + t65 * (t86 * t40 + t47 * t41);
t90 = t51 / 0.2e1;
t88 = pkin(3) * t47;
t87 = pkin(3) * t48;
t5 = t33 ^ 2 - t54 ^ 2;
t83 = t5 * qJD(2);
t46 = -t50 * pkin(3) - pkin(2);
t78 = qJD(2) * t46;
t49 = sin(qJ(2));
t77 = qJD(2) * t49;
t76 = qJD(2) * t50;
t75 = qJD(2) * t51;
t74 = qJD(3) * t48;
t73 = qJD(3) * t50;
t72 = qJD(3) * t51;
t71 = qJD(4) * t46;
t14 = t33 * t87 + t46 * t54;
t70 = t14 * qJD(2);
t15 = -t46 * t33 + t54 * t87;
t69 = t15 * qJD(2);
t42 = -t48 ^ 2 + t50 ^ 2;
t66 = t42 * qJD(2);
t64 = pkin(2) * t48 * qJD(2);
t63 = pkin(2) * t76;
t62 = t33 * t78;
t61 = t54 * t78;
t60 = t48 * t76;
t57 = t86 * qJD(3);
t56 = t86 * qJD(4);
t55 = t49 * t65;
t21 = t65 * t54;
t22 = t54 * t33 * qJD(2);
t20 = t65 * t33;
t19 = t53 * t51 - t54 * t90;
t18 = t33 * t90 + t52 * t51;
t10 = t17 * qJD(2);
t8 = t16 * qJD(2);
t2 = t19 * qJD(2) + t33 * t55;
t1 = t18 * qJD(2) + t54 * t55;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t77, -t75, 0, 0, 0, 0, 0, -t48 * t72 - t49 * t76, t48 * t77 - t50 * t72, 0, 0, 0, 0, 0, t65 * t19 + t33 * t77, t65 * t18 + t54 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t75 - t49 * t73, t49 * t74 - t50 * t75, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t16, -t65 * t17; 0, 0, 0, 0, t48 * t73, t42 * qJD(3), 0, 0, 0, -pkin(2) * t74, -pkin(2) * t73, -t33 * t21, t65 * t5, 0, 0, 0, t14 * qJD(3) + t54 * t71, t15 * qJD(3) - t33 * t71; 0, 0, 0, 0, t60, t66, t73, -t74, 0, -pkin(5) * t73 - t64, pkin(5) * t74 - t63, -t22, t83, -t20, -t21, 0, t70 + t92, t69 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t83, -t20, -t21, 0, t61 + t92, -t62 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10; 0, 0, 0, 0, -t60, -t66, 0, 0, 0, t64, t63, t22, -t83, 0, 0, 0, t68 - t70, t67 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t88, -pkin(3) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t88, (-t57 - t56) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t83, 0, 0, 0, t68 - t61, t67 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t88, pkin(3) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
