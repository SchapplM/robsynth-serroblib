% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPPR5
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
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% tau_reg [4x6]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:11
% EndTime: 2018-11-14 14:05:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (43->20), mult. (75->30), div. (0->0), fcn. (68->4), ass. (0->13)
t11 = qJD(4) ^ 2;
t7 = sin(pkin(5));
t15 = t7 ^ 2 * qJDD(1) - g(1);
t14 = t7 * qJDD(1);
t13 = qJDD(1) - g(1);
t10 = cos(qJ(4));
t8 = cos(pkin(5));
t9 = sin(qJ(4));
t12 = t8 * t10 + t7 * t9;
t2 = t7 * t10 - t8 * t9;
t6 = qJDD(2) + g(3);
t3 = -t8 * qJDD(1) + qJDD(3);
t1 = [t13, t8 ^ 2 * qJDD(1) + t15, -t3 * t8 + t15, 0, -qJDD(4) * t12 - t2 * t11, -t2 * qJDD(4) + t12 * t11; 0, t6, t6, 0, 0, 0; 0, 0, -g(2) * t7 - t13 * t8 + qJDD(3), 0, t10 * qJDD(4) - t11 * t9, -qJDD(4) * t9 - t11 * t10; 0, 0, 0, qJDD(4), g(1) * t12 - g(2) * t2 + t10 * t3 - t9 * t14, g(1) * t2 + g(2) * t12 - t10 * t14 - t9 * t3;];
tau_reg  = t1;
