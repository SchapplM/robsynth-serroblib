% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:49
% EndTime: 2019-12-05 16:00:51
% DurationCPUTime: 0.55s
% Computational Cost: add. (354->76), mult. (785->110), div. (0->0), fcn. (790->6), ass. (0->65)
t62 = qJD(4) + qJD(5);
t53 = -pkin(2) - pkin(6);
t93 = -pkin(7) + t53;
t48 = sin(qJ(4));
t50 = cos(qJ(5));
t86 = t50 * t48;
t47 = sin(qJ(5));
t51 = cos(qJ(4));
t87 = t47 * t51;
t34 = t86 + t87;
t23 = t62 * t34;
t36 = t93 * t48;
t37 = t93 * t51;
t85 = t50 * t51;
t88 = t47 * t48;
t33 = -t85 + t88;
t49 = sin(qJ(2));
t55 = -t88 / 0.2e1 + t85 / 0.2e1;
t16 = (t33 / 0.2e1 + t55) * t49;
t74 = t16 * qJD(1);
t92 = -t74 + t62 * (-t50 * t36 - t47 * t37);
t54 = -t86 / 0.2e1 - t87 / 0.2e1;
t17 = (t34 / 0.2e1 + t54) * t49;
t73 = t17 * qJD(1);
t91 = -t73 + t62 * (t47 * t36 - t50 * t37);
t90 = -t49 / 0.2e1;
t89 = pkin(4) * t51;
t84 = pkin(4) * qJD(5);
t83 = qJD(4) * pkin(4);
t7 = -t33 ^ 2 + t34 ^ 2;
t78 = t7 * qJD(2);
t43 = pkin(4) * t48 + qJ(3);
t77 = qJD(5) * t43;
t14 = -t33 * t43 + t34 * t89;
t76 = t14 * qJD(2);
t15 = -t33 * t89 - t34 * t43;
t75 = t15 * qJD(2);
t72 = t33 * qJD(2);
t71 = t34 * qJD(2);
t42 = t48 ^ 2 - t51 ^ 2;
t70 = t42 * qJD(2);
t69 = t48 * qJD(2);
t68 = t48 * qJD(4);
t45 = t49 * qJD(2);
t67 = t51 * qJD(2);
t66 = t51 * qJD(4);
t52 = cos(qJ(2));
t65 = t52 * qJD(2);
t64 = qJ(3) * qJD(4);
t63 = qJD(2) * qJ(3);
t61 = t43 * t72;
t60 = t43 * t71;
t59 = t48 * t67;
t58 = t48 * t63;
t57 = t51 * t63;
t56 = pkin(4) * t62;
t22 = t62 * t33;
t24 = t33 * t71;
t19 = t33 * t90 + t55 * t49;
t18 = t34 * t90 + t54 * t49;
t10 = t17 * qJD(2);
t8 = t16 * qJD(2);
t2 = t19 * qJD(2) + t52 * t23;
t1 = t18 * qJD(2) - t52 * t22;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t45, -t65, t45, t65, (-pkin(2) * t49 + qJ(3) * t52) * qJD(2) + t49 * qJD(3), 0, 0, 0, 0, 0, t48 * t65 + t49 * t66, -t49 * t68 + t51 * t65, 0, 0, 0, 0, 0, t62 * t19 + t34 * t65, t62 * t18 - t33 * t65; 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t67 + t52 * t68, -t48 * t45 + t52 * t66, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t16, -t62 * t17; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t48 * t66, t42 * qJD(4), 0, 0, 0, qJD(3) * t48 + t51 * t64, qJD(3) * t51 - t48 * t64, t34 * t22, t62 * t7, 0, 0, 0, qJD(3) * t34 + qJD(4) * t14 - t33 * t77, -qJD(3) * t33 + qJD(4) * t15 - t34 * t77; 0, 0, 0, 0, 0, qJD(2), t63, 0, 0, 0, 0, 0, t69, t67, 0, 0, 0, 0, 0, t71, -t72; 0, 0, 0, 0, 0, 0, 0, -t59, t70, -t68, -t66, 0, -t53 * t68 + t57, -t53 * t66 - t58, t24, t78, -t23, t22, 0, t76 + t92, t75 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t78, -t23, t22, 0, -t61 + t92, -t60 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t63, 0, 0, 0, 0, 0, -t69, -t67, 0, 0, 0, 0, 0, -t71, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t66, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10; 0, 0, 0, 0, 0, 0, 0, t59, -t70, 0, 0, 0, -t57, t58, -t24, -t78, 0, 0, 0, t74 - t76, t73 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t84, -t50 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t56, -t50 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t78, 0, 0, 0, t74 + t61, t73 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t83, t50 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
