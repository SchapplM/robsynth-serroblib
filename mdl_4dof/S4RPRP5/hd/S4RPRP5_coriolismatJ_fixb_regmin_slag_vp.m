% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:03
% EndTime: 2019-12-31 16:45:04
% DurationCPUTime: 0.26s
% Computational Cost: add. (395->57), mult. (799->78), div. (0->0), fcn. (870->4), ass. (0->47)
t38 = sin(pkin(6));
t39 = cos(pkin(6));
t40 = sin(qJ(3));
t61 = cos(qJ(3));
t29 = t61 * t38 + t40 * t39;
t25 = t29 ^ 2;
t60 = pkin(5) + qJ(2);
t33 = t38 ^ 2 + t39 ^ 2;
t27 = t40 * t38 - t61 * t39;
t10 = pkin(3) * t29 + t27 * qJ(4);
t35 = -t39 * pkin(2) - pkin(1);
t41 = t27 * pkin(3) - t29 * qJ(4);
t8 = t35 + t41;
t1 = t8 * t10;
t59 = t1 * qJD(1);
t3 = t10 * t27 + t8 * t29;
t58 = t3 * qJD(1);
t4 = -t10 * t29 + t8 * t27;
t57 = t4 * qJD(1);
t32 = t60 * t39;
t42 = t60 * t38;
t13 = t40 * t32 + t61 * t42;
t14 = t61 * t32 - t40 * t42;
t5 = t13 * t29 - t14 * t27;
t56 = t5 * qJD(1);
t55 = t10 * qJD(1);
t24 = t27 ^ 2;
t7 = t24 - t25;
t54 = t7 * qJD(1);
t12 = t24 + t25;
t53 = t12 * qJD(1);
t52 = t13 * qJD(3);
t11 = t14 * qJD(3);
t51 = t25 * qJD(1);
t50 = t27 * qJD(1);
t18 = t27 * qJD(3);
t49 = t29 * qJD(1);
t20 = t29 * qJD(3);
t48 = t29 * qJD(4);
t31 = t33 * qJ(2);
t47 = t31 * qJD(1);
t46 = t33 * qJD(1);
t45 = qJD(3) * qJ(4);
t44 = t27 * t49;
t43 = t35 * t49;
t22 = t29 * qJD(2);
t2 = [0, 0, 0, 0, 0, t33 * qJD(2), t31 * qJD(2), -t27 * t20, t7 * qJD(3), 0, 0, 0, t35 * t20, -t35 * t18, t3 * qJD(3) - t27 * t48, t12 * qJD(2), t4 * qJD(3) + t25 * qJD(4), t5 * qJD(2) + t1 * qJD(3) - t8 * t48; 0, 0, 0, 0, 0, t46, t47, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t56; 0, 0, 0, 0, 0, 0, 0, -t44, t54, -t18, -t20, 0, -t11 + t43, -t35 * t50 + t52, -t11 + t58, qJD(3) * t41 - t27 * qJD(4), -t52 + t57, t59 + (-t14 * pkin(3) - t13 * qJ(4)) * qJD(3) + t14 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t18, t51, -t8 * t49 + t11; 0, 0, 0, 0, 0, -t46, -t47, 0, 0, 0, 0, 0, t20, -t18, t20, -t53, t18, qJD(3) * t10 - t48 - t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t50, t49, 0, t50, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, t44, -t54, 0, 0, 0, -t22 - t43, (qJD(1) * t35 + qJD(2)) * t27, -t22 - t58, 0, -t27 * qJD(2) - t57, -qJD(2) * t10 - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t50, -t49, 0, -t50, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t51, (qJD(1) * t8 + qJD(2)) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
