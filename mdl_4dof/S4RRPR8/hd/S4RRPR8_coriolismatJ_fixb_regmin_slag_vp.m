% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:25
% EndTime: 2019-12-31 17:08:26
% DurationCPUTime: 0.38s
% Computational Cost: add. (300->70), mult. (584->97), div. (0->0), fcn. (546->4), ass. (0->68)
t55 = qJD(2) - qJD(4);
t37 = sin(qJ(4));
t38 = sin(qJ(2));
t39 = cos(qJ(4));
t40 = cos(qJ(2));
t18 = t38 * t37 + t40 * t39;
t20 = -t40 * t37 + t38 * t39;
t4 = t18 ^ 2 - t20 ^ 2;
t88 = t4 * qJD(1);
t87 = pkin(5) - pkin(6);
t86 = t55 * t18;
t43 = t55 * t20;
t28 = t87 * t38;
t29 = t87 * t40;
t85 = t55 * (t37 * t28 + t39 * t29);
t84 = t55 * (t39 * t28 - t37 * t29);
t83 = pkin(2) + pkin(3);
t78 = t38 * qJ(3);
t11 = t83 * t40 + pkin(1) + t78;
t82 = t11 * t18;
t76 = t40 * qJ(3);
t16 = -t83 * t38 + t76;
t1 = -t11 * t20 + t16 * t18;
t81 = t1 * qJD(1);
t2 = t16 * t20 + t82;
t80 = t2 * qJD(1);
t42 = -t40 * pkin(2) - t78;
t24 = -pkin(1) + t42;
t27 = t38 * pkin(2) - t76;
t7 = t24 * t40 + t27 * t38;
t71 = t7 * qJD(1);
t8 = -t24 * t38 + t27 * t40;
t70 = t8 * qJD(1);
t69 = qJD(1) * t20;
t68 = qJD(1) * t38;
t67 = qJD(1) * t40;
t66 = qJD(2) * t38;
t34 = qJD(2) * t40;
t65 = qJD(3) * t38;
t64 = qJD(4) * t11;
t36 = t38 ^ 2;
t30 = t40 ^ 2 - t36;
t63 = t30 * qJD(1);
t62 = t36 * qJD(1);
t61 = t37 * qJD(2);
t60 = t37 * qJD(3);
t59 = t39 * qJD(2);
t58 = t39 * qJD(3);
t57 = t40 * qJD(3);
t56 = qJD(2) * qJ(3);
t54 = pkin(1) * t68;
t53 = pkin(1) * t67;
t52 = pkin(5) * t66;
t51 = pkin(5) * t34;
t50 = qJD(1) * t82;
t49 = t11 * t69;
t48 = t18 * t69;
t47 = t24 * t27 * qJD(1);
t46 = t24 * t68;
t45 = t18 * t68;
t44 = t20 * t68;
t41 = t42 * qJD(2) + t57;
t31 = t38 * t67;
t26 = t55 * t39;
t25 = t55 * t37;
t23 = t39 * qJ(3) - t37 * t83;
t22 = t37 * qJ(3) + t39 * t83;
t3 = [0, 0, 0, t38 * t34, t30 * qJD(2), 0, 0, 0, -pkin(1) * t66, -pkin(1) * t34, -t8 * qJD(2) + t38 * t57, 0, -t7 * qJD(2) + t36 * qJD(3), (qJD(2) * t27 - t65) * t24, t18 * t43, -t55 * t4, 0, 0, 0, t1 * qJD(2) + t18 * t65 + t20 * t64, t2 * qJD(2) - t18 * t64 + t20 * t65; 0, 0, 0, t31, t63, t34, -t66, 0, -t51 - t54, t52 - t53, -t51 - t70, t41, -t52 - t71, t41 * pkin(5) + t47, t48, -t88, -t86, -t43, 0, t81 - t85, t80 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34, t62, -t46 + t51, 0, 0, 0, 0, 0, t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t88, t86, t43, 0, t49 + t85, -t50 + t84; 0, 0, 0, -t31, -t63, 0, 0, 0, t54, t53, t70, 0, t71, -t47, -t48, t88, 0, 0, 0, -t81, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, 0, 0, 0, 0, t23 * qJD(4) + t60, -t22 * qJD(4) + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t56, 0, 0, 0, 0, 0, t61, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t23, -t55 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t62, t46, 0, 0, 0, 0, 0, -t45, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t56, 0, 0, 0, 0, 0, -t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t88, 0, 0, 0, -t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23 * qJD(2) - t60, t22 * qJD(2) - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
