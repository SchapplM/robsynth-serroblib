% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S1R1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
%
% Output:
% m_new_reg [(3*2)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S1R1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_invdynm_fixb_reg2_snew_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_invdynm_fixb_reg2_snew_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1R1_invdynm_fixb_reg2_snew_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1R1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_invdynm_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:57
% EndTime: 2020-06-19 09:12:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (36->19), mult. (72->19), div. (0->0), fcn. (56->2), ass. (0->12)
t25 = sin(qJ(1));
t26 = cos(qJ(1));
t22 = t25 * g(1) - t26 * g(2);
t23 = t26 * g(1) + t25 * g(2);
t30 = -t25 * t22 - t26 * t23;
t27 = qJD(1) ^ 2;
t21 = t26 * qJDD(1) - t25 * t27;
t29 = -pkin(1) * t21 - t25 * g(3);
t28 = t26 * t22 - t25 * t23;
t20 = t25 * qJDD(1) + t26 * t27;
t18 = -pkin(1) * t20 + t26 * g(3);
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t21, 0, -t20, 0, t29, -t18, -t28, -pkin(1) * t28; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t20, 0, t21, 0, t18, t29, t30, pkin(1) * t30; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t22, t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t27, 0, 0, -g(3), -t22, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, qJDD(1), 0, g(3), 0, -t23, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t22, t23, 0, 0;];
m_new_reg = t1;
