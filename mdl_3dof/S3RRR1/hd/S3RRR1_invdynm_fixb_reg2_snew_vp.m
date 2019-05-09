% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3RRR1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3RRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:32:49
% EndTime: 2019-05-04 18:32:50
% DurationCPUTime: 0.96s
% Computational Cost: add. (3023->120), mult. (4196->132), div. (0->0), fcn. (2694->6), ass. (0->78)
t246 = qJD(1) + qJD(2);
t242 = qJD(3) + t246;
t240 = t242 ^ 2;
t245 = qJDD(1) + qJDD(2);
t241 = qJDD(3) + t245;
t247 = sin(qJ(3));
t250 = cos(qJ(3));
t223 = t250 * t240 + t247 * t241;
t226 = t247 * t240 - t250 * t241;
t248 = sin(qJ(2));
t251 = cos(qJ(2));
t197 = t251 * t223 - t248 * t226;
t214 = pkin(5) * t223 - t250 * g(3);
t280 = pkin(5) * t226 - t247 * g(3);
t183 = pkin(4) * t197 + t251 * t214 - t248 * t280;
t249 = sin(qJ(1));
t252 = cos(qJ(1));
t201 = t248 * t223 + t251 * t226;
t279 = t252 * t197 - t249 * t201;
t290 = pkin(4) * t201 + t248 * t214 + t251 * t280;
t302 = pkin(3) * t279 + t252 * t183 - t249 * t290;
t292 = t249 * t197 + t252 * t201;
t301 = pkin(3) * t292 + t249 * t183 + t252 * t290;
t238 = t249 * g(1) - t252 * g(2);
t234 = qJDD(1) * pkin(1) + t238;
t239 = t252 * g(1) + t249 * g(2);
t254 = qJD(1) ^ 2;
t235 = -t254 * pkin(1) - t239;
t209 = -t251 * t234 + t248 * t235;
t204 = t245 * pkin(2) - t209;
t210 = t248 * t234 + t251 * t235;
t244 = t246 ^ 2;
t205 = -t244 * pkin(2) + t210;
t187 = -t250 * t204 + t247 * t205;
t188 = t247 * t204 + t250 * t205;
t263 = t247 * t187 + t250 * t188;
t176 = t250 * t187 - t247 * t188;
t265 = t251 * t176;
t170 = -t248 * t263 + t265;
t267 = t248 * t176;
t284 = t251 * t263 + t267;
t299 = t249 * t170 + t252 * t284;
t298 = t252 * t170 - t249 * t284;
t229 = t251 * t244 + t248 * t245;
t218 = pkin(4) * t229 - t251 * g(3);
t232 = t248 * t244 - t251 * t245;
t257 = t252 * t229 - t249 * t232;
t281 = pkin(4) * t232 - t248 * g(3);
t293 = pkin(3) * t257 + t252 * t218 - t249 * t281;
t278 = t249 * t229 + t252 * t232;
t291 = pkin(3) * t278 + t249 * t218 + t252 * t281;
t262 = t248 * t209 + t251 * t210;
t193 = t251 * t209 - t248 * t210;
t264 = t252 * t193;
t286 = -t249 * t262 + t264;
t266 = t249 * t193;
t285 = t252 * t262 + t266;
t260 = -t249 * t238 - t252 * t239;
t237 = t252 * qJDD(1) - t249 * t254;
t259 = -pkin(3) * t237 - t249 * g(3);
t258 = -pkin(2) * t226 - t187;
t256 = t252 * t238 - t249 * t239;
t255 = -pkin(2) * t223 - t188;
t253 = pkin(1) * g(3);
t236 = t249 * qJDD(1) + t252 * t254;
t227 = -pkin(3) * t236 + t252 * g(3);
t196 = -pkin(1) * t229 - t210;
t195 = -pkin(1) * t232 - t209;
t190 = pkin(1) * t193;
t189 = pkin(4) * t262 + t253;
t179 = -pkin(1) * t197 + t255;
t178 = -pkin(1) * t201 + t258;
t173 = pkin(2) * t176;
t172 = pkin(2) * g(3) + pkin(5) * t263;
t167 = -pkin(1) * t170 - t173;
t166 = pkin(4) * t170 + pkin(5) * t265 - t248 * t172;
t165 = pkin(4) * t284 + pkin(5) * t267 + t251 * t172 + t253;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t237, 0, -t236, 0, t259, -t227, -t256, -pkin(3) * t256, 0, 0, -t278, 0, -t257, 0, t291, t293, t286, pkin(3) * t286 + pkin(4) * t264 - t249 * t189, 0, 0, -t292, 0, -t279, 0, t301, t302, t298, pkin(3) * t298 - t249 * t165 + t252 * t166; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t236, 0, t237, 0, t227, t259, t260, pkin(3) * t260, 0, 0, t257, 0, -t278, 0, -t293, t291, t285, pkin(3) * t285 + pkin(4) * t266 + t252 * t189, 0, 0, t279, 0, -t292, 0, -t302, t301, t299, pkin(3) * t299 + t252 * t165 + t249 * t166; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t238, t239, 0, 0, 0, 0, 0, 0, 0, t245, t195, t196, 0, -t190, 0, 0, 0, 0, 0, t241, t178, t179, 0, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t254, 0, 0, -g(3), -t238, 0, 0, 0, -t232, 0, -t229, 0, t281, t218, t193, pkin(4) * t193, 0, 0, -t201, 0, -t197, 0, t290, t183, t170, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, 0, qJDD(1), 0, g(3), 0, -t239, 0, 0, 0, t229, 0, -t232, 0, -t218, t281, t262, t189, 0, 0, t197, 0, -t201, 0, -t183, t290, t284, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t238, t239, 0, 0, 0, 0, 0, 0, 0, t245, t195, t196, 0, -t190, 0, 0, 0, 0, 0, t241, t178, t179, 0, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, 0, -t244, 0, 0, -g(3), t209, 0, 0, 0, -t226, 0, -t223, 0, t280, t214, t176, pkin(5) * t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, t245, 0, g(3), 0, t210, 0, 0, 0, t223, 0, -t226, 0, -t214, t280, t263, t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, -t209, -t210, 0, 0, 0, 0, 0, 0, 0, t241, t258, t255, 0, -t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, -t240, 0, 0, -g(3), t187, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, t241, 0, g(3), 0, t188, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t187, -t188, 0, 0;];
m_new_reg  = t1;
