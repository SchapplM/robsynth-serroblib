% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% tau_reg [4x6]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:15
% EndTime: 2018-11-14 13:38:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (31->16), mult. (48->21), div. (0->0), fcn. (42->4), ass. (0->12)
t11 = qJD(4) ^ 2;
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t8 = cos(pkin(5));
t7 = sin(pkin(5));
t6 = qJDD(1) - g(3);
t5 = -t9 * qJDD(4) - t10 * t11;
t4 = -t10 * qJDD(4) + t9 * t11;
t3 = -g(1) * t7 + g(2) * t8 + qJDD(2);
t2 = t7 * t10 + t8 * t9;
t1 = t8 * t10 - t7 * t9;
t12 = [t6, t6, t6, 0, 0, 0; 0, t3, t3, 0, t5, t4; 0, 0, -g(1) * t8 - g(2) * t7 + qJDD(3), 0, -t4, t5; 0, 0, 0, qJDD(4), -g(1) * t1 - g(2) * t2 - t9 * qJDD(2) + t10 * qJDD(3), g(1) * t2 - g(2) * t1 - t10 * qJDD(2) - t9 * qJDD(3);];
tau_reg  = t12;
