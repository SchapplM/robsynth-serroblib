% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% tau_reg [2x6]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S2RR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (36->18), mult. (48->23), div. (0->0), fcn. (26->6), ass. (0->14)
t9 = qJ(1) + qJ(2);
t5 = sin(t9);
t6 = cos(t9);
t17 = g(1) * t6 + g(2) * t5;
t12 = cos(qJ(2));
t16 = t12 * qJDD(1) * pkin(1) + g(1) * t5 - g(2) * t6;
t8 = qJD(1) + qJD(2);
t15 = qJD(1) * (-qJD(2) + t8);
t14 = qJD(2) * (-qJD(1) - t8);
t13 = cos(qJ(1));
t11 = sin(qJ(1));
t10 = sin(qJ(2));
t7 = qJDD(1) + qJDD(2);
t1 = [qJDD(1), g(1) * t11 - g(2) * t13, g(1) * t13 + g(2) * t11, t7, (t10 * t14 + t12 * t7) * pkin(1) + t16, ((-qJDD(1) - t7) * t10 + t12 * t14) * pkin(1) + t17; 0, 0, 0, t7, t10 * pkin(1) * t15 + t16, (-qJDD(1) * t10 + t12 * t15) * pkin(1) + t17;];
tau_reg = t1;
