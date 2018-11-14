% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S3PPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% tau_reg [3x5]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3PPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:16
% EndTime: 2018-11-14 10:11:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (21->7), mult. (34->9), div. (0->0), fcn. (24->2), ass. (0->8)
t3 = qJDD(1) - g(1);
t7 = qJDD(2) + g(2);
t6 = qJD(3) ^ 2;
t5 = cos(qJ(3));
t4 = sin(qJ(3));
t2 = -t4 * qJDD(3) - t5 * t6;
t1 = -t5 * qJDD(3) + t4 * t6;
t8 = [t3, t3, 0, t2, t1; 0, t7, 0, -t1, t2; 0, 0, qJDD(3), -t3 * t4 + t7 * t5, -t3 * t5 - t7 * t4;];
tau_reg  = t8;
