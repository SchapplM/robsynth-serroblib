% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:43
% EndTime: 2019-12-05 15:28:45
% DurationCPUTime: 0.38s
% Computational Cost: add. (629->60), mult. (1307->78), div. (0->0), fcn. (1430->4), ass. (0->48)
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t58 = sin(qJ(4));
t83 = cos(qJ(4));
t45 = t58 * t56 - t83 * t57;
t47 = t83 * t56 + t58 * t57;
t20 = pkin(4) * t47 + t45 * qJ(5);
t66 = t47 * qJD(5);
t84 = -t20 * qJD(4) + t66;
t42 = t45 ^ 2;
t43 = t47 ^ 2;
t82 = pkin(6) + qJ(3);
t51 = t56 ^ 2 + t57 ^ 2;
t53 = -t57 * pkin(3) - pkin(2);
t59 = t45 * pkin(4) - t47 * qJ(5);
t18 = t53 + t59;
t2 = t18 * t20;
t81 = t2 * qJD(2);
t5 = t18 * t47 + t20 * t45;
t80 = t5 * qJD(2);
t6 = t18 * t45 - t20 * t47;
t79 = t6 * qJD(2);
t50 = t82 * t57;
t60 = t82 * t56;
t27 = t58 * t50 + t83 * t60;
t28 = t83 * t50 - t58 * t60;
t9 = t27 * t47 - t45 * t28;
t76 = t9 * qJD(2);
t75 = t20 * qJD(2);
t17 = t42 - t43;
t74 = t17 * qJD(2);
t22 = t42 + t43;
t73 = t22 * qJD(2);
t72 = t27 * qJD(4);
t21 = t28 * qJD(4);
t71 = t43 * qJD(2);
t70 = t45 * qJD(2);
t69 = t45 * qJD(4);
t68 = t47 * qJD(2);
t67 = t47 * qJD(4);
t49 = t51 * qJ(3);
t65 = t49 * qJD(2);
t64 = t51 * qJD(2);
t63 = qJD(4) * qJ(5);
t62 = t45 * t68;
t61 = t53 * t68;
t40 = t47 * qJD(3);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t69, -t67, 0, -t69, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t51 * qJD(3), t49 * qJD(3), -t45 * t67, t17 * qJD(4), 0, 0, 0, t53 * t67, -t53 * t69, t5 * qJD(4) - t45 * t66, t22 * qJD(3), t6 * qJD(4) + t43 * qJD(5), t9 * qJD(3) + t2 * qJD(4) - t18 * t66; 0, 0, 0, 0, 0, 0, t64, t65, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, -t62, t74, -t69, -t67, 0, -t21 + t61, -t53 * t70 + t72, -t21 + t80, t59 * qJD(4) - t45 * qJD(5), -t72 + t79, t81 + (-t28 * pkin(4) - t27 * qJ(5)) * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t69, t71, -t18 * t68 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t64, -t65, 0, 0, 0, 0, 0, t67, -t69, t67, -t73, t69, -t76 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t70, t68, 0, t70, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t62, -t74, 0, 0, 0, -t40 - t61, (qJD(2) * t53 + qJD(3)) * t45, -t40 - t80, 0, -t45 * qJD(3) - t79, -qJD(3) * t20 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t70, -t68, 0, -t70, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t71, (qJD(2) * t18 + qJD(3)) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
