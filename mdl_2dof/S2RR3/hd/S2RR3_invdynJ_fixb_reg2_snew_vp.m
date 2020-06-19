% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [2x(2*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S2RR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (43->14), mult. (69->21), div. (0->0), fcn. (46->4), ass. (0->14)
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t15 = t11 * g(1) - t13 * g(2);
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t4 = qJDD(1) * pkin(1) + t15;
t14 = t13 * g(1) + t11 * g(2);
t5 = -qJD(1) ^ 2 * pkin(1) - t14;
t1 = -t10 * t5 + t12 * t4;
t2 = t10 * t4 + t12 * t5;
t9 = qJD(1) + qJD(2);
t8 = qJDD(1) + qJDD(2);
t7 = t9 ^ 2;
t3 = [0, 0, 0, 0, 0, qJDD(1), t15, t14, 0, 0, 0, 0, 0, 0, 0, t8, pkin(1) * (-t10 * t7 + t12 * t8) + t1, pkin(1) * (-t10 * t8 - t12 * t7) - t2, 0, pkin(1) * (t12 * t1 + t10 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t1, -t2, 0, 0;];
tauJ_reg = t3;
