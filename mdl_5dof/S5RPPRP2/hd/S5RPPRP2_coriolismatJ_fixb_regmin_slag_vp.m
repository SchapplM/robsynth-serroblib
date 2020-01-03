% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:26
% EndTime: 2019-12-31 17:49:27
% DurationCPUTime: 0.41s
% Computational Cost: add. (777->62), mult. (1455->79), div. (0->0), fcn. (1578->6), ass. (0->50)
t58 = cos(pkin(8));
t57 = sin(pkin(8));
t60 = sin(qJ(4));
t84 = t60 * t57;
t85 = cos(qJ(4));
t47 = -t85 * t58 + t84;
t62 = t85 * t57;
t49 = t60 * t58 + t62;
t23 = pkin(4) * t49 + t47 * qJ(5);
t67 = t49 * qJD(5);
t87 = -qJD(4) * t23 + t67;
t45 = t47 ^ 2;
t46 = t49 ^ 2;
t53 = sin(pkin(7)) * pkin(1) + qJ(3);
t86 = pkin(6) + t53;
t52 = t57 ^ 2 + t58 ^ 2;
t51 = -cos(pkin(7)) * pkin(1) - t58 * pkin(3) - pkin(2);
t61 = t47 * pkin(4) - t49 * qJ(5);
t18 = t51 + t61;
t2 = t18 * t23;
t83 = t2 * qJD(1);
t33 = t86 * t58;
t20 = t60 * t33 + t86 * t62;
t21 = t85 * t33 - t86 * t84;
t5 = t20 * t49 - t21 * t47;
t82 = t5 * qJD(1);
t6 = t18 * t49 + t23 * t47;
t81 = t6 * qJD(1);
t7 = t18 * t47 - t23 * t49;
t80 = t7 * qJD(1);
t77 = t23 * qJD(1);
t17 = t45 - t46;
t76 = t17 * qJD(1);
t75 = t20 * qJD(4);
t19 = t21 * qJD(4);
t24 = t45 + t46;
t74 = t24 * qJD(1);
t29 = t52 * t53;
t73 = t29 * qJD(1);
t72 = t46 * qJD(1);
t71 = t47 * qJD(1);
t70 = t47 * qJD(4);
t69 = t49 * qJD(1);
t68 = t49 * qJD(4);
t66 = t52 * qJD(1);
t65 = qJD(4) * qJ(5);
t64 = t47 * t69;
t63 = t51 * t69;
t43 = t49 * qJD(3);
t1 = [0, 0, 0, 0, 0, 0, t52 * qJD(3), t29 * qJD(3), -t47 * t68, t17 * qJD(4), 0, 0, 0, t51 * t68, -t51 * t70, t6 * qJD(4) - t47 * t67, t24 * qJD(3), t7 * qJD(4) + t46 * qJD(5), t5 * qJD(3) + t2 * qJD(4) - t18 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t66, t73, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, -t64, t76, -t70, -t68, 0, -t19 + t63, -t51 * t71 + t75, -t19 + t81, qJD(4) * t61 - t47 * qJD(5), -t75 + t80, t83 + (-t21 * pkin(4) - t20 * qJ(5)) * qJD(4) + t21 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t70, t72, -t18 * t69 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t70, -t68, 0, -t70, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, -t66, -t73, 0, 0, 0, 0, 0, t68, -t70, t68, -t74, t70, -t82 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t71, t69, 0, t71, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, t64, -t76, 0, 0, 0, -t43 - t63, (qJD(1) * t51 + qJD(3)) * t47, -t43 - t81, 0, -t47 * qJD(3) - t80, -qJD(3) * t23 - t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t71, -t69, 0, -t71, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t72, (qJD(1) * t18 + qJD(3)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
