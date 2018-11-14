% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:30
% EndTime: 2018-11-14 13:53:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (293->52), mult. (674->82), div. (0->0), fcn. (546->6), ass. (0->52)
t38 = sin(pkin(7));
t64 = t38 * pkin(2);
t41 = sin(qJ(2));
t63 = t38 * t41;
t39 = cos(pkin(7));
t62 = t39 * t41;
t43 = cos(qJ(2));
t29 = (-t38 * t43 - t62) * pkin(1);
t40 = sin(qJ(4));
t61 = t40 * t29;
t30 = (t39 * t43 - t63) * pkin(1);
t60 = t40 * t30;
t42 = cos(qJ(4));
t59 = t42 * t30;
t58 = pkin(1) * qJD(1);
t57 = pkin(1) * qJD(2);
t37 = t43 * pkin(1) + pkin(2);
t27 = pkin(1) * t62 + t38 * t37;
t49 = -pkin(1) * t63 + t39 * t37;
t5 = t27 * t30 + t49 * t29;
t56 = t5 * qJD(1);
t25 = pkin(3) + t49;
t8 = -t42 * t25 + t40 * t27;
t55 = t8 * qJD(1);
t9 = t40 * t25 + t42 * t27;
t54 = t9 * qJD(1);
t24 = t42 * t29;
t12 = -t24 + t60;
t53 = t12 * qJD(1);
t13 = t59 + t61;
t52 = t13 * qJD(1);
t36 = t39 * pkin(2) + pkin(3);
t51 = -t36 / 0.2e1 - t25 / 0.2e1;
t50 = pkin(1) * (-qJD(1) - qJD(2));
t48 = t64 / 0.2e1 + t27 / 0.2e1;
t47 = t30 / 0.2e1 + t51;
t1 = t47 * t42 + (t29 / 0.2e1 + t48) * t40;
t26 = -t42 * t36 + t40 * t64;
t46 = -t1 * qJD(1) - t26 * qJD(2);
t28 = t40 * t36 + t42 * t64;
t44 = t48 * t42;
t3 = -t24 / 0.2e1 - t44 + t47 * t40;
t45 = -t3 * qJD(1) + t28 * qJD(2);
t23 = t28 * qJD(4);
t22 = t26 * qJD(4);
t11 = t13 * qJD(2);
t10 = t12 * qJD(2);
t7 = t9 * qJD(4);
t6 = t8 * qJD(4);
t4 = -t60 / 0.2e1 + t24 / 0.2e1 - t44 + t51 * t40;
t2 = -t59 / 0.2e1 - t61 / 0.2e1 + t51 * t42 + t48 * t40;
t14 = [0, 0, 0, 0, -t41 * t57, -t43 * t57, t5 * qJD(2), 0, -t10 - t7, -t11 + t6; 0, 0, 0, 0, t41 * t50, t43 * t50, t56 + (t29 * t39 + t30 * t38) * qJD(2) * pkin(2), 0, t4 * qJD(4) - t10 - t53, t2 * qJD(4) - t11 - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t4 * qJD(2) - t54 - t7, t2 * qJD(2) + t55 + t6; 0, 0, 0, 0, t41 * t58, t43 * t58, -t56, 0, t3 * qJD(4) + t53, t1 * qJD(4) + t52; 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t23 - t45, t22 - t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t3 * qJD(2) + t54, -t1 * qJD(2) - t55; 0, 0, 0, 0, 0, 0, 0, 0, t45, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t14;
