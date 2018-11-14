% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:38
% EndTime: 2018-11-14 13:43:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (192->57), mult. (188->65), div. (0->0), fcn. (82->6), ass. (0->40)
t27 = pkin(6) + qJ(2);
t25 = qJ(3) + t27;
t15 = sin(t25);
t16 = cos(t25);
t49 = -g(1) * t16 - g(2) * t15;
t48 = t16 * pkin(3) + t15 * qJ(4);
t26 = qJDD(2) + qJDD(3);
t47 = t26 * pkin(3);
t30 = sin(qJ(3));
t31 = cos(qJ(3));
t44 = pkin(2) * qJD(2);
t36 = qJD(3) * t44;
t43 = pkin(2) * qJDD(2);
t46 = t30 * t43 + t31 * t36;
t45 = -t30 * t36 + t31 * t43;
t42 = qJD(3) * t30;
t41 = qJD(3) * t31;
t40 = t30 * t44;
t39 = t31 * t44;
t28 = qJD(2) + qJD(3);
t38 = t28 * t42;
t37 = -t15 * pkin(3) + t16 * qJ(4);
t35 = t46 + t49;
t21 = t26 * qJ(4);
t22 = t28 * qJD(4);
t1 = t21 + t22 + t46;
t34 = g(1) * t15 - g(2) * t16 + t45;
t33 = -qJDD(4) + t34;
t32 = t28 * t39 - t35;
t29 = qJDD(1) - g(3);
t24 = cos(t27);
t23 = sin(t27);
t18 = -t31 * pkin(2) - pkin(3);
t17 = t30 * pkin(2) + qJ(4);
t8 = pkin(2) * t41 + qJD(4);
t5 = t28 * t40;
t4 = t28 * qJ(4) + t40;
t3 = -t28 * pkin(3) + qJD(4) - t39;
t2 = qJDD(4) - t45 - t47;
t6 = [t29, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, qJDD(2), g(1) * t23 - g(2) * t24, g(1) * t24 + g(2) * t23, t26 (t26 * t31 - t38) * pkin(2) + t34 (-t26 * t30 - t28 * t41) * pkin(2) - t35, -pkin(2) * t38 + (pkin(3) - t18) * t26 + t33, t17 * t26 + t8 * t28 + t1 + t49, t1 * t17 + t4 * t8 + t2 * t18 + t3 * pkin(2) * t42 - g(1) * (-pkin(2) * t23 + t37) - g(2) * (pkin(2) * t24 + t48); 0, 0, 0, 0, t26, t34 + t5, t32, t33 + t5 + 0.2e1 * t47, 0.2e1 * t21 + 0.2e1 * t22 - t32, t1 * qJ(4) + t4 * qJD(4) - t2 * pkin(3) - g(1) * t37 - g(2) * t48 + (-t3 * t30 - t31 * t4) * t44; 0, 0, 0, 0, 0, 0, 0, -t26, -t28 ^ 2, -t4 * t28 - t33 - t47;];
tau_reg  = t6;
