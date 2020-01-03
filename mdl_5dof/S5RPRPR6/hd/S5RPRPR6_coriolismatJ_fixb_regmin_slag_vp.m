% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:51
% EndTime: 2019-12-31 18:17:52
% DurationCPUTime: 0.34s
% Computational Cost: add. (299->74), mult. (564->73), div. (0->0), fcn. (440->6), ass. (0->54)
t42 = sin(qJ(3));
t49 = cos(pkin(8)) * pkin(1) + pkin(2);
t70 = cos(qJ(3));
t72 = pkin(1) * sin(pkin(8));
t26 = t42 * t49 + t70 * t72;
t24 = qJ(4) + t26;
t77 = -qJ(4) / 0.2e1 - t24 / 0.2e1;
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t33 = t41 ^ 2 - t43 ^ 2;
t38 = qJD(1) + qJD(3);
t76 = t38 * t33;
t60 = t26 * qJD(3);
t61 = t26 * qJD(1);
t75 = t61 + t60;
t36 = qJD(4) * t41;
t25 = t42 * t72 - t70 * t49;
t62 = t25 * qJD(3);
t69 = -t41 * t62 + t36;
t37 = qJD(4) * t43;
t68 = -t43 * t62 + t37;
t47 = -pkin(3) + t25;
t1 = -t24 * t25 + t47 * t26;
t67 = t1 * qJD(1);
t66 = qJD(5) * t41;
t65 = qJD(5) * t43;
t64 = qJD(5) * (-pkin(3) - pkin(7));
t63 = t25 * qJD(1);
t59 = t41 * qJD(1);
t58 = t43 * qJD(1);
t57 = -t62 + qJD(4);
t56 = qJ(4) * qJD(3);
t55 = qJ(4) * qJD(5);
t54 = t24 * t59;
t53 = t24 * t58;
t52 = t25 * t59;
t51 = t25 * t58;
t50 = t41 * t65;
t29 = t38 * t43;
t48 = t26 / 0.2e1 + t77;
t2 = t48 * t41;
t46 = -t2 * qJD(1) + t41 * t56;
t3 = t48 * t43;
t45 = t3 * qJD(1) - t43 * t56;
t39 = qJ(4) * qJD(4);
t32 = t38 * qJ(4);
t31 = t33 * qJD(5);
t28 = t38 * t41;
t27 = t41 * t29;
t23 = -pkin(7) + t47;
t12 = t24 * qJD(4);
t5 = t24 * t43;
t4 = (-t26 / 0.2e1 + t77) * t41;
t6 = [0, 0, 0, 0, 0, -t60, t62, t60, t57, t1 * qJD(3) + t12, -t50, t31, 0, 0, 0, t24 * t65 + t69, -t24 * t66 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t75, t63 + t62, t75, -t63 + t57, t67 + (-t26 * pkin(3) - t25 * qJ(4)) * qJD(3) + t12, -t50, t31, 0, 0, 0, t5 * qJD(5) - t52 + t69, t4 * qJD(5) - t51 + t68; 0, 0, 0, 0, 0, 0, 0, 0, t38, t38 * t24, 0, 0, 0, 0, 0, t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t76, -t66, -t65, 0, t5 * qJD(3) - t23 * t66 + t53, t4 * qJD(3) - t23 * t65 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t66; 0, 0, 0, 0, 0, t61, -t63, -t61, t63 + qJD(4), t39 - t67, -t50, t31, 0, 0, 0, -t3 * qJD(5) + t36 + t52, t2 * qJD(5) + t37 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t39, -t50, t31, 0, 0, 0, t43 * t55 + t36, -t41 * t55 + t37; 0, 0, 0, 0, 0, 0, 0, 0, t38, t32, 0, 0, 0, 0, 0, t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t76, -t66, -t65, 0, -t41 * t64 - t45, -t43 * t64 - t46; 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t24 * qJD(1) - t56, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t32, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t76, 0, 0, 0, t3 * qJD(3) - t53, -t2 * qJD(3) + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t76, 0, 0, 0, t45, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
