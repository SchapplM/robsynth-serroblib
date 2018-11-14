% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:40
% EndTime: 2018-11-14 13:51:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (124->43), mult. (324->61), div. (0->0), fcn. (168->4), ass. (0->36)
t22 = sin(pkin(6));
t25 = cos(qJ(2));
t36 = pkin(1) * qJD(2);
t23 = cos(pkin(6));
t24 = sin(qJ(2));
t40 = t23 * t24;
t10 = (t22 * t25 + t40) * t36;
t42 = qJD(1) * t10;
t41 = t22 * t24;
t39 = t23 * t25;
t21 = qJD(1) + qJD(2);
t37 = pkin(1) * qJD(1);
t33 = t25 * t37;
t13 = t21 * pkin(2) + t33;
t34 = t24 * t37;
t16 = t23 * t34;
t5 = t22 * t13 + t16;
t19 = t25 * pkin(1) + pkin(2);
t38 = pkin(1) * t40 + t22 * t19;
t35 = pkin(1) * t41;
t32 = t22 * t34;
t17 = t36 * t39;
t31 = (-qJD(2) + t21) * t37;
t30 = (-qJD(1) - t21) * t36;
t28 = -qJD(2) * t35 + t17;
t27 = t23 * t19 - t35;
t7 = qJD(1) * t17 - qJD(2) * t32;
t4 = t23 * t13 - t32;
t20 = t21 * qJD(4);
t3 = t20 + t7;
t11 = (t39 - t41) * t37;
t9 = t22 * t33 + t16;
t8 = qJD(4) + t28;
t2 = t21 * qJ(4) + t5;
t1 = -t21 * pkin(3) + qJD(4) - t4;
t6 = [0, 0, 0, 0, t24 * t30, t25 * t30, -t4 * t10 - t27 * t42 + t5 * t28 + t7 * t38, -t10 * t21 - t42, t8 * t21 + t3, t3 * (qJ(4) + t38) + t2 * t8 + t42 * (-pkin(3) - t27) + t1 * t10; 0, 0, 0, 0, t24 * t31, t25 * t31, -t5 * t11 + t4 * t9 + (t22 * t7 - t23 * t42) * pkin(2), t9 * t21 - t42, -t11 * t21 + 0.2e1 * t20 + t7, t3 * (t22 * pkin(2) + qJ(4)) + t42 * (-t23 * pkin(2) - pkin(3)) - t1 * t9 + (qJD(4) - t11) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t21 ^ 2, -t2 * t21 + t42;];
tauc_reg  = t6;
