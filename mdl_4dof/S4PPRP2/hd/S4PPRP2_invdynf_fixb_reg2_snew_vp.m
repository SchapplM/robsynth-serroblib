% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRP2
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
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:44:37
% EndTime: 2019-05-04 18:44:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (437->61), mult. (614->43), div. (0->0), fcn. (472->4), ass. (0->28)
t306 = sin(qJ(3));
t307 = cos(qJ(3));
t308 = qJD(3) ^ 2;
t297 = t306 * qJDD(3) + t307 * t308;
t298 = t307 * qJDD(3) - t306 * t308;
t304 = sin(pkin(5));
t305 = cos(pkin(5));
t310 = -t304 * t297 + t305 * t298;
t284 = t305 * t297 + t304 * t298;
t302 = -g(2) + qJDD(1);
t295 = t304 * g(1) + t305 * t302;
t296 = -t305 * g(1) + t304 * t302;
t283 = t306 * t295 + t307 * t296;
t282 = t307 * t295 - t306 * t296;
t301 = g(3) - qJDD(2);
t281 = -t304 * t295 + t305 * t296;
t280 = t305 * t295 + t304 * t296;
t279 = -qJDD(3) * pkin(3) - t308 * qJ(4) + qJDD(4) - t282;
t278 = -t308 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t283;
t277 = -t306 * t282 + t307 * t283;
t276 = t307 * t282 + t306 * t283;
t275 = t307 * t278 + t306 * t279;
t274 = t306 * t278 - t307 * t279;
t273 = -t304 * t276 + t305 * t277;
t272 = t305 * t276 + t304 * t277;
t271 = -t304 * t274 + t305 * t275;
t270 = t305 * t274 + t304 * t275;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, 0, 0, 0, 0, 0, 0, -t284, -t310, 0, t273, 0, 0, 0, 0, 0, 0, -t284, 0, t310, t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, 0, 0, 0, 0, 0, 0, t310, -t284, 0, t272, 0, 0, 0, 0, 0, 0, t310, 0, t284, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, 0, 0, 0, 0, 0, 0, -t284, -t310, 0, t273, 0, 0, 0, 0, 0, 0, -t284, 0, t310, t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, 0, 0, 0, 0, 0, 0, t310, -t284, 0, t272, 0, 0, 0, 0, 0, 0, t310, 0, t284, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, 0, 0, 0, 0, 0, 0, -t297, -t298, 0, t277, 0, 0, 0, 0, 0, 0, -t297, 0, t298, t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, 0, 0, 0, 0, 0, 0, t298, -t297, 0, t276, 0, 0, 0, 0, 0, 0, t298, 0, t297, t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, -qJDD(3), 0, t283, 0, 0, 0, 0, 0, 0, -t308, 0, qJDD(3), t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t308, 0, t282, 0, 0, 0, 0, 0, 0, qJDD(3), 0, t308, -t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, 0, qJDD(3), t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3), 0, -t308, t279;];
f_new_reg  = t1;
