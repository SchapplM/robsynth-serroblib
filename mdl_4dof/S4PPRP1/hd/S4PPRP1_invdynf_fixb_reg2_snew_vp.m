% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:43:17
% EndTime: 2019-05-04 18:43:17
% DurationCPUTime: 0.37s
% Computational Cost: add. (318->67), mult. (502->45), div. (0->0), fcn. (384->4), ass. (0->24)
t258 = sin(qJ(3));
t259 = qJD(3) ^ 2;
t262 = cos(qJ(3));
t248 = t258 * qJDD(3) + t262 * t259;
t249 = -t262 * qJDD(3) + t258 * t259;
t256 = sin(pkin(5));
t257 = cos(pkin(5));
t263 = t257 * t248 + t256 * t249;
t236 = -t256 * t248 + t257 * t249;
t250 = t256 * g(1) - t257 * g(2);
t247 = -qJDD(2) + t250;
t251 = -t257 * g(1) - t256 * g(2);
t239 = -t258 * t247 + t262 * t251;
t238 = -t262 * t247 - t258 * t251;
t254 = g(3) - qJDD(1);
t245 = t257 * t251;
t244 = t256 * t251;
t235 = -qJDD(3) * pkin(3) - t259 * qJ(4) + qJDD(4) - t238;
t234 = -t259 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t239;
t233 = -t258 * t238 + t262 * t239;
t232 = t262 * t238 + t258 * t239;
t231 = t262 * t234 + t258 * t235;
t230 = t258 * t234 - t262 * t235;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256 * t250 + t245, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256 * t247 + t245, 0, 0, 0, 0, 0, 0, -t263, t236, 0, t256 * t232 + t257 * t233, 0, 0, 0, 0, 0, 0, -t263, 0, -t236, t256 * t230 + t257 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t257 * t250 + t244, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257 * t247 + t244, 0, 0, 0, 0, 0, 0, t236, t263, 0, -t257 * t232 + t256 * t233, 0, 0, 0, 0, 0, 0, t236, 0, -t263, -t257 * t230 + t256 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, -t248, t249, 0, t233, 0, 0, 0, 0, 0, 0, -t248, 0, -t249, t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, 0, 0, 0, 0, 0, 0, t249, t248, 0, -t232, 0, 0, 0, 0, 0, 0, t249, 0, -t248, -t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, -t248, t249, 0, t233, 0, 0, 0, 0, 0, 0, -t248, 0, -t249, t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, 0, 0, 0, 0, 0, 0, -t249, -t248, 0, t232, 0, 0, 0, 0, 0, 0, -t249, 0, t248, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, -qJDD(3), 0, t239, 0, 0, 0, 0, 0, 0, -t259, 0, qJDD(3), t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t259, 0, t238, 0, 0, 0, 0, 0, 0, qJDD(3), 0, t259, -t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, 0, qJDD(3), t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3), 0, -t259, t235;];
f_new_reg  = t1;
