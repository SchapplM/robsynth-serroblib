% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:07:14
% EndTime: 2018-11-14 14:07:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (125->48), mult. (236->56), div. (0->0), fcn. (184->6), ass. (0->36)
t26 = sin(qJ(3));
t27 = cos(qJ(3));
t37 = qJD(1) * qJD(3);
t46 = qJDD(1) * t26 + t27 * t37;
t22 = pkin(5) + qJ(3);
t20 = sin(t22);
t21 = cos(t22);
t45 = g(1) * t21 - g(2) * t20;
t25 = cos(pkin(5));
t24 = sin(pkin(5));
t42 = t26 * t24;
t10 = -t27 * t25 + t42;
t6 = t10 * qJD(1);
t44 = qJD(4) + t6;
t41 = t26 * t25;
t40 = qJDD(3) * pkin(3);
t38 = qJDD(1) * t27;
t36 = qJDD(3) * qJ(4);
t35 = t24 * t38 + t46 * t25;
t34 = qJD(1) * t42;
t32 = t46 * t24 - t25 * t38 + t37 * t41;
t31 = -t6 + t34;
t11 = t27 * t24 + t41;
t9 = t11 * qJD(3);
t30 = -t9 * qJD(3) - t10 * qJDD(3);
t7 = t11 * qJD(1);
t29 = g(1) * t20 + g(2) * t21 - t35;
t2 = qJDD(4) + t32 - t40;
t28 = t7 * qJD(3) - t32 - t45;
t23 = qJDD(2) + g(3);
t8 = t10 * qJD(3);
t5 = qJD(3) * qJ(4) + t7;
t4 = -qJD(3) * pkin(3) + t44;
t3 = -t8 * qJD(3) + t11 * qJDD(3);
t1 = t36 + (qJD(4) - t34) * qJD(3) + t35;
t12 = [qJDD(1) - g(1), -g(1) + (t24 ^ 2 + t25 ^ 2) * qJDD(1), 0, t30, -t3, t30, t3, t1 * t11 + t2 * t10 + t4 * t9 - t5 * t8 - g(1); 0, t23, 0, 0, 0, 0, 0, t23; 0, 0, qJDD(3), t28, t31 * qJD(3) + t29, -qJDD(4) + t28 + 0.2e1 * t40, 0.2e1 * t36 + (0.2e1 * qJD(4) - t31) * qJD(3) - t29, t1 * qJ(4) - t2 * pkin(3) - t4 * t7 - g(1) * (t21 * pkin(3) + t20 * qJ(4)) - g(2) * (-t20 * pkin(3) + t21 * qJ(4)) + t44 * t5; 0, 0, 0, 0, 0, -qJDD(3), -qJD(3) ^ 2, -t5 * qJD(3) + t2 + t45;];
tau_reg  = t12;
