% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPPR2
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
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:38:09
% EndTime: 2019-05-04 18:38:10
% DurationCPUTime: 0.28s
% Computational Cost: add. (285->44), mult. (338->31), div. (0->0), fcn. (312->4), ass. (0->27)
t210 = sin(qJ(4));
t211 = cos(qJ(4));
t212 = qJD(4) ^ 2;
t204 = t210 * qJDD(4) + t211 * t212;
t205 = -t211 * qJDD(4) + t210 * t212;
t208 = sin(pkin(5));
t209 = cos(pkin(5));
t213 = -t208 * t204 + t209 * t205;
t207 = -g(2) + qJDD(1);
t202 = t208 * g(1) + t209 * t207;
t194 = t209 * t204 + t208 * t205;
t206 = g(3) - qJDD(2);
t203 = -t209 * g(1) + t208 * t207;
t201 = -qJDD(3) + t202;
t199 = t209 * t203;
t198 = t208 * t203;
t193 = -t208 * t202 + t199;
t192 = t209 * t202 + t198;
t191 = -t210 * t201 + t211 * t203;
t190 = -t211 * t201 - t210 * t203;
t189 = -t208 * t201 + t199;
t188 = t209 * t201 + t198;
t187 = -t210 * t190 + t211 * t191;
t186 = t211 * t190 + t210 * t191;
t185 = t208 * t186 + t209 * t187;
t184 = -t209 * t186 + t208 * t187;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, 0, 0, 0, 0, 0, -t194, t213, 0, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, 0, 0, 0, 0, 0, t213, t194, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, 0, 0, 0, 0, 0, -t194, t213, 0, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, 0, 0, 0, 0, 0, t213, t194, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, -t204, t205, 0, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, 0, 0, 0, 0, 0, 0, t205, t204, 0, -t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, -t204, t205, 0, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, 0, 0, 0, 0, 0, 0, -t205, -t204, 0, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, -qJDD(4), 0, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t212, 0, t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206;];
f_new_reg  = t1;
