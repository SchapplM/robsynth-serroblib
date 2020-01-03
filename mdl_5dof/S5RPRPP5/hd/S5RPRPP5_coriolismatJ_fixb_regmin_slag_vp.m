% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:44
% EndTime: 2019-12-31 18:16:46
% DurationCPUTime: 0.50s
% Computational Cost: add. (355->101), mult. (551->100), div. (0->0), fcn. (396->2), ass. (0->71)
t47 = cos(qJ(3));
t41 = t47 * qJ(4);
t46 = sin(qJ(3));
t83 = pkin(3) + pkin(4);
t85 = -t83 * t46 + t41;
t24 = -t46 * pkin(3) + t41;
t33 = t46 * qJD(4);
t52 = t24 * qJD(3) + t33;
t84 = t83 / 0.2e1;
t81 = t47 * pkin(3);
t80 = t85 * qJD(3) + t33;
t14 = qJ(2) - t85;
t76 = t46 * qJ(4);
t16 = -t83 * t47 - t76;
t1 = t14 * t16;
t79 = t1 * qJD(1);
t2 = -t14 * t47 + t16 * t46;
t78 = t2 * qJD(1);
t3 = t14 * t46 + t16 * t47;
t77 = t3 * qJD(1);
t49 = -pkin(1) - pkin(6);
t72 = qJ(5) + t49;
t18 = t72 * t46;
t19 = t72 * t47;
t6 = t18 * t46 + t19 * t47;
t75 = t6 * qJD(1);
t20 = qJ(2) - t24;
t23 = t76 + t81;
t8 = t20 * t47 + t23 * t46;
t74 = t8 * qJD(1);
t9 = -t20 * t46 + t23 * t47;
t73 = t9 * qJD(1);
t71 = qJD(2) * t47;
t10 = t76 + (t84 + pkin(3) / 0.2e1 + pkin(4) / 0.2e1) * t47;
t70 = t10 * qJD(1);
t69 = t14 * qJD(1);
t68 = t18 * qJD(3);
t67 = t20 * qJD(1);
t44 = t46 ^ 2;
t45 = t47 ^ 2;
t56 = t44 / 0.2e1 + t45 / 0.2e1;
t21 = 0.1e1 / 0.2e1 + t56;
t66 = t21 * qJD(1);
t25 = t44 + t45;
t65 = t25 * qJD(1);
t26 = t44 - t45;
t64 = t26 * qJD(1);
t34 = t46 * qJD(3);
t39 = t47 * qJD(1);
t38 = t47 * qJD(3);
t63 = t47 * qJD(4);
t62 = qJ(2) * qJD(3);
t61 = qJD(1) * qJ(2);
t60 = t23 * t67;
t59 = t20 * t39;
t58 = t49 * t34;
t57 = t49 * t38;
t55 = t46 * t61;
t54 = t47 * t61;
t53 = t45 * qJD(4) - t71;
t51 = qJD(2) - t63;
t40 = qJD(2) * t46;
t50 = -t46 * t63 + t40;
t43 = qJ(4) * qJD(4);
t42 = qJD(3) * qJ(4);
t37 = t45 * qJD(1);
t35 = t46 * qJD(1);
t27 = t46 * t39;
t22 = -0.1e1 / 0.2e1 + t56;
t11 = -t81 / 0.2e1 + (t84 - pkin(4) / 0.2e1) * t47;
t4 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t46 * t38, t26 * qJD(3), 0, 0, 0, t47 * t62 + t40, -t46 * t62 + t71, t8 * qJD(3) + t50, 0, -t9 * qJD(3) + t53, (qJD(3) * t23 + t51) * t20, -t2 * qJD(3) + t50, t3 * qJD(3) + t53, t25 * qJD(5), -t1 * qJD(3) + t6 * qJD(5) + t51 * t14; 0, 0, 0, 0, qJD(1), t61, 0, 0, 0, 0, 0, t35, t39, t35, 0, -t39, t67, t35, -t39, 0, t22 * qJD(5) + t69; 0, 0, 0, 0, 0, 0, -t27, t64, -t34, -t38, 0, t54 - t58, -t55 - t57, -t58 + t74, -t52, t57 - t73, t52 * t49 + t60, -t68 - t78, t19 * qJD(3) + t77, t80, -t79 + (t19 * qJ(4) - t18 * t83) * qJD(3) + t18 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t34, t37, t58 - t59, -t27, t37, t34, -t14 * t39 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t22 * qJD(2) + t11 * qJD(3) + t75; 0, 0, 0, 0, -qJD(1), -t61, 0, 0, 0, 0, 0, -t35, -t39, -t35, 0, t39, -t67, -t35, t39, 0, t21 * qJD(5) - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t38, -t34, 0, t38, t52, -t34, t38, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66; 0, 0, 0, 0, 0, 0, t27, -t64, 0, 0, 0, -t54, t55, -t74, 0, t73, -t60, t47 * qJD(5) + t78, t46 * qJD(5) - t77, 0, t10 * qJD(5) + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t43, 0, qJD(4), 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t42, 0, qJD(3), 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t35, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t37, t59, t27, -t37, 0, (-qJD(5) + t69) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t42, 0, -qJD(3), 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t34, -t65, -t21 * qJD(2) - t10 * qJD(3) + t63 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t35, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
