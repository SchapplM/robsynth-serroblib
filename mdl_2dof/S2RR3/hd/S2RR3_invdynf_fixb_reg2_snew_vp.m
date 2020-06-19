% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*3)x(3*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S2RR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynf_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynf_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynf_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:29
% EndTime: 2020-06-19 09:14:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (147->37), mult. (218->36), div. (0->0), fcn. (156->4), ass. (0->23)
t152 = qJD(1) + qJD(2);
t150 = t152 ^ 2;
t151 = qJDD(1) + qJDD(2);
t153 = sin(qJ(2));
t155 = cos(qJ(2));
t141 = t153 * t150 - t155 * t151;
t154 = sin(qJ(1));
t156 = cos(qJ(1));
t158 = -t155 * t150 - t153 * t151;
t162 = t154 * t141 + t156 * t158;
t161 = t156 * t141 - t154 * t158;
t147 = t154 * g(1) - t156 * g(2);
t148 = -t156 * g(1) - t154 * g(2);
t157 = qJD(1) ^ 2;
t146 = -t154 * qJDD(1) - t156 * t157;
t145 = t156 * qJDD(1) - t154 * t157;
t144 = -t157 * pkin(1) + t148;
t143 = qJDD(1) * pkin(1) + t147;
t138 = t153 * t143 + t155 * t144;
t137 = t155 * t143 - t153 * t144;
t136 = -t153 * t137 + t155 * t138;
t135 = t155 * t137 + t153 * t138;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t146, -t145, 0, -t154 * t147 + t156 * t148, 0, 0, 0, 0, 0, 0, t162, t161, 0, -t154 * t135 + t156 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t145, t146, 0, t156 * t147 + t154 * t148, 0, 0, 0, 0, 0, 0, -t161, t162, 0, t156 * t135 + t154 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, -qJDD(1), 0, t148, 0, 0, 0, 0, 0, 0, t158, t141, 0, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t157, 0, t147, 0, 0, 0, 0, 0, 0, -t141, t158, 0, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t151, 0, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t150, 0, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
f_new_reg = t1;
