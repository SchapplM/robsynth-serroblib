% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:35:39
% EndTime: 2019-05-04 18:35:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (210->56), mult. (286->33), div. (0->0), fcn. (264->4), ass. (0->23)
t200 = sin(pkin(5));
t201 = cos(pkin(5));
t195 = t200 * g(1) - t201 * g(2);
t191 = -qJDD(2) + t195;
t208 = t200 * t191;
t202 = sin(qJ(4));
t203 = cos(qJ(4));
t204 = qJD(4) ^ 2;
t193 = t203 * qJDD(4) - t202 * t204;
t194 = t202 * qJDD(4) + t203 * t204;
t207 = t201 * t193 - t200 * t194;
t196 = t201 * g(1) + t200 * g(2);
t192 = -qJDD(3) + t196;
t206 = t202 * t191 - t203 * t192;
t205 = t200 * t193 + t201 * t194;
t199 = g(3) - qJDD(1);
t190 = t201 * t196;
t189 = t200 * t196;
t187 = t201 * t191;
t186 = -t203 * t191 - t202 * t192;
t184 = t203 * t186 - t202 * t206;
t183 = t202 * t186 + t203 * t206;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200 * t195 - t190, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190 - t208, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201 * t192 - t208, 0, 0, 0, 0, 0, 0, t207, -t205, 0, t201 * t183 + t200 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t201 * t195 - t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189 + t187, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200 * t192 + t187, 0, 0, 0, 0, 0, 0, t205, t207, 0, t200 * t183 - t201 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, 0, 0, 0, 0, 0, 0, t193, -t194, 0, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, t194, t193, 0, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, -t193, t194, 0, -t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, 0, 0, 0, 0, 0, 0, -t194, -t193, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, 0, 0, 0, 0, 0, 0, -t194, -t193, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, 0, 0, 0, 0, 0, 0, t193, -t194, 0, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, -qJDD(4), 0, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t204, 0, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199;];
f_new_reg  = t1;
