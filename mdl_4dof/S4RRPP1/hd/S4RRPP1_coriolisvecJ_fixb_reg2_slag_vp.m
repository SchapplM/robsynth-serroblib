% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:40
% EndTime: 2018-11-14 13:51:41
% DurationCPUTime: 0.17s
% Computational Cost: add. (140->43), mult. (384->64), div. (0->0), fcn. (200->4), ass. (0->40)
t23 = sin(pkin(6));
t25 = sin(qJ(2));
t44 = t23 * t25;
t24 = cos(pkin(6));
t43 = t24 * t25;
t26 = cos(qJ(2));
t42 = t24 * t26;
t39 = pkin(1) * qJD(2);
t18 = t39 * t42;
t40 = pkin(1) * qJD(1);
t37 = t25 * t40;
t36 = t23 * t37;
t7 = qJD(1) * t18 - qJD(2) * t36;
t20 = t26 * pkin(1) + pkin(2);
t41 = pkin(1) * t43 + t23 * t20;
t38 = pkin(1) * t44;
t22 = qJD(1) + qJD(2);
t21 = t22 * qJD(4);
t3 = t21 + t7;
t11 = (t42 - t44) * t40;
t35 = t11 * t22 - t7;
t34 = (-qJD(2) + t22) * t40;
t33 = (-qJD(1) - t22) * t39;
t32 = t23 * t26 + t43;
t14 = t22 * pkin(2) + t26 * t40;
t5 = t23 * t14 + t24 * t37;
t12 = -qJD(2) * t38 + t18;
t31 = t24 * t20 - t38;
t30 = pkin(1) * t32;
t4 = t24 * t14 - t36;
t10 = qJD(2) * t30;
t29 = t32 * qJD(1) * t39;
t6 = qJD(1) * t10;
t9 = qJD(1) * t30;
t28 = t9 * t22 - t29;
t27 = -t10 * t22 - t29;
t8 = qJD(4) + t12;
t2 = t22 * qJ(4) + t5;
t1 = -t22 * pkin(3) + qJD(4) - t4;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t33, t26 * t33, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t12 * t22 - t7, 0, -t4 * t10 + t5 * t12 - t6 * t31 + t7 * t41, 0, 0, 0, 0, 0, 0, t27, 0, t8 * t22 + t3, t3 * (qJ(4) + t41) + t2 * t8 + t6 * (-pkin(3) - t31) + t1 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t34, t26 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t28, t35, 0, -t5 * t11 + t4 * t9 + (t23 * t7 - t6 * t24) * pkin(2), 0, 0, 0, 0, 0, 0, t28, 0, 0.2e1 * t21 - t35, t3 * (t23 * pkin(2) + qJ(4)) + t6 * (-t24 * pkin(2) - pkin(3)) - t1 * t9 + (qJD(4) - t11) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 ^ 2, -t2 * t22 + t6;];
tauc_reg  = t13;
