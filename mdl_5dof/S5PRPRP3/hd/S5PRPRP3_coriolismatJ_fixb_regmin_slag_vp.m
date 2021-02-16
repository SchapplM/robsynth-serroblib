% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP3
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
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:03
% EndTime: 2021-01-15 15:14:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (420->57), mult. (938->91), div. (0->0), fcn. (926->6), ass. (0->57)
t60 = sin(qJ(4));
t58 = t60 ^ 2;
t61 = cos(qJ(4));
t59 = t61 ^ 2;
t50 = t58 + t59;
t75 = t50 * qJD(2);
t87 = t61 * pkin(4);
t86 = cos(qJ(2));
t85 = sin(qJ(2));
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t45 = t79 * t85 - t80 * t86;
t26 = t45 * t60;
t56 = -t80 * pkin(2) - pkin(3);
t47 = t56 - t87;
t84 = t47 * t60;
t52 = t79 * t86;
t53 = t80 * t85;
t46 = -t52 - t53;
t3 = (-0.1e1 + t50) * t46 * t45;
t83 = t3 * qJD(1);
t4 = pkin(4) * t84;
t82 = t4 * qJD(2);
t55 = t79 * pkin(2) + pkin(6);
t81 = qJ(5) + t55;
t30 = t60 * t87 - t84;
t78 = t30 * qJD(2);
t35 = t58 * pkin(4) + t47 * t61;
t77 = t35 * qJD(2);
t43 = t81 * t61;
t76 = t43 * qJD(4);
t51 = t59 - t58;
t74 = t51 * qJD(2);
t73 = t60 * qJD(2);
t72 = t60 * qJD(4);
t71 = t60 * qJD(5);
t70 = t61 * qJD(2);
t57 = t61 * qJD(4);
t69 = pkin(4) * t72;
t68 = pkin(4) * t73;
t67 = t56 * t73;
t66 = t56 * t70;
t65 = t60 * t70;
t64 = t46 * t57;
t63 = t53 / 0.2e1 + t52 / 0.2e1;
t42 = t81 * t60;
t23 = t42 * t60 + t43 * t61;
t10 = (t59 / 0.2e1 + t58 / 0.2e1) * t46 + t63;
t62 = -t10 * qJD(1) + t23 * qJD(2);
t27 = t45 * t61;
t11 = t63 - t50 * t46 / 0.2e1;
t8 = t27 * qJD(2) - t46 * t72;
t7 = t26 * qJD(2) + t64;
t6 = t26 * qJD(4) + t46 * t70;
t5 = t27 * qJD(4) - t46 * t73;
t2 = pkin(4) * t26;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(2); 0, 0, -t85 * qJD(2), -t86 * qJD(2), (-t79 * t45 + t80 * t46) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, t6, t5, t6, t5, -t45 * t75, t83 + (-t23 * t45 - t46 * t47) * qJD(2) + t2 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, t8, 0, pkin(4) * t64 + t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * qJD(5) - t83; 0, 0, 0, 0, 0, t60 * t57, t51 * qJD(4), 0, 0, 0, t56 * t72, t56 * t57, -t30 * qJD(4), t35 * qJD(4), t50 * qJD(5), t4 * qJD(4) + t23 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t65, t74, t57, -t72, 0, -t55 * t57 + t67, t55 * t72 + t66, -t76 - t78, t42 * qJD(4) + t77, -pkin(4) * t57, -pkin(4) * t76 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t57, -t72, -t57, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t65, -t74, 0, 0, 0, -t67, -t66, -t71 + t78, -t61 * qJD(5) - t77, 0, -pkin(4) * t71 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t70, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t57, -t75, -t62 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t70, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
