% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRP4
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% tau_reg [4x6]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:14
% EndTime: 2018-11-14 14:06:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (63->29), mult. (105->34), div. (0->0), fcn. (70->2), ass. (0->15)
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t18 = -t11 * qJD(1) + t12 * qJD(2);
t17 = g(2) * t12;
t10 = qJDD(1) - g(1);
t15 = qJDD(2) + g(2);
t7 = t12 * qJD(1) + t11 * qJD(2);
t13 = qJD(3) ^ 2;
t8 = t12 * qJDD(2);
t5 = -t11 * qJDD(3) - t12 * t13;
t4 = -t12 * qJDD(3) + t11 * t13;
t3 = qJD(3) * pkin(3) + t18;
t2 = qJD(3) * t18 + t12 * qJDD(1) + t11 * qJDD(2);
t1 = qJDD(3) * pkin(3) - t7 * qJD(3) - t11 * qJDD(1) + t8;
t6 = [t10, t10, 0, t5, t4, -t1 * t11 + t2 * t12 - g(1) + (-t11 * t7 - t12 * t3) * qJD(3); 0, t15, 0, -t4, t5, t1 * t12 + t2 * t11 + g(2) + (-t11 * t3 + t12 * t7) * qJD(3); 0, 0, qJDD(3), -t10 * t11 + t17 + t8, -t10 * t12 - t15 * t11 (t3 - t18) * t7 + (g(1) * t11 + t1 + t17) * pkin(3); 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t6;
