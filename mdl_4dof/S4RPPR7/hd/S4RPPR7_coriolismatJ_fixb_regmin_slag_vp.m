% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:42
% EndTime: 2019-12-31 16:41:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (98->41), mult. (209->52), div. (0->0), fcn. (204->4), ass. (0->35)
t27 = -pkin(1) - qJ(3);
t42 = -pkin(5) + t27;
t41 = sin(qJ(4));
t25 = sin(pkin(6));
t21 = t25 ^ 2;
t26 = cos(pkin(6));
t22 = t26 ^ 2;
t16 = t21 + t22;
t28 = cos(qJ(4));
t11 = -t41 * t25 + t28 * t26;
t9 = t28 * t25 + t41 * t26;
t1 = -t11 ^ 2 + t9 ^ 2;
t40 = t1 * qJD(1);
t8 = t16 * t27;
t39 = t8 * qJD(1);
t38 = t9 * qJD(1);
t3 = t9 * qJD(4);
t17 = t25 * pkin(3) + qJ(2);
t37 = qJD(1) * t17;
t36 = t11 * qJD(1);
t6 = t11 * qJD(4);
t30 = -t21 / 0.2e1 - t22 / 0.2e1;
t15 = -0.1e1 / 0.2e1 + t30;
t35 = t15 * qJD(1);
t34 = t16 * qJD(1);
t33 = t25 * qJD(1);
t32 = t26 * qJD(1);
t31 = t9 * t36;
t29 = qJD(3) + t37;
t24 = qJ(2) * qJD(2);
t23 = qJD(1) * qJ(2);
t14 = 0.1e1 / 0.2e1 + t30;
t13 = t42 * t26;
t12 = t42 * t25;
t2 = [0, 0, 0, 0, qJD(2), t24, qJD(2) * t25, qJD(2) * t26, t16 * qJD(3), -t8 * qJD(3) + t24, -t9 * t6, t1 * qJD(4), 0, 0, 0, qJD(2) * t9 + t17 * t6, qJD(2) * t11 - t17 * t3; 0, 0, 0, 0, qJD(1), t23, t33, t32, 0, t14 * qJD(3) + t23, 0, 0, 0, 0, 0, t38, t36; 0, 0, 0, 0, 0, 0, 0, 0, t34, t14 * qJD(2) - t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t40, -t3, -t6, 0, t17 * t36 + (-t28 * t12 - t41 * t13) * qJD(4), -t9 * t37 + (t41 * t12 - t28 * t13) * qJD(4); 0, 0, 0, 0, -qJD(1), -t23, -t33, -t32, 0, t15 * qJD(3) - t23, 0, 0, 0, 0, 0, -t38, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t6; 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t15 * qJD(2) + t39, 0, 0, 0, 0, 0, t6, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t40, 0, 0, 0, -t29 * t11, t29 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
