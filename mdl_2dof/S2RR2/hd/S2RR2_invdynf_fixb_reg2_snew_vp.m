% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S2RR2
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% f_new_reg [(3*3)x(3*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S2RR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_invdynf_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_invdynf_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_invdynf_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:19:16
% EndTime: 2019-05-04 18:19:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (127->45), mult. (306->62), div. (0->0), fcn. (194->4), ass. (0->36)
t201 = sin(qJ(2));
t199 = t201 ^ 2;
t203 = cos(qJ(2));
t200 = t203 ^ 2;
t210 = t199 + t200;
t209 = qJD(1) * qJD(2);
t208 = t201 * qJDD(1);
t207 = t203 * qJDD(1);
t202 = sin(qJ(1));
t204 = cos(qJ(1));
t193 = -t204 * g(1) + t202 * g(3);
t192 = t202 * g(1) + t204 * g(3);
t206 = qJD(1) ^ 2;
t205 = qJD(2) ^ 2;
t196 = t203 * t206 * t201;
t195 = -t200 * t206 - t205;
t194 = -t199 * t206 - t205;
t191 = -qJDD(2) + t196;
t190 = qJDD(2) + t196;
t189 = t210 * t206;
t188 = -t204 * qJDD(1) + t202 * t206;
t187 = t202 * qJDD(1) + t204 * t206;
t186 = t210 * qJDD(1);
t185 = -0.2e1 * t201 * t209 + t207;
t184 = 0.2e1 * t203 * t209 + t208;
t183 = t206 * pkin(1) + t192;
t182 = qJDD(1) * pkin(1) + t193;
t181 = -t201 * g(2) + t203 * t182;
t180 = -t203 * g(2) - t201 * t182;
t179 = t203 * t191 - t201 * t194;
t178 = -t201 * t190 + t203 * t195;
t177 = t201 * t191 + t203 * t194;
t176 = t203 * t190 + t201 * t195;
t175 = -t201 * t180 + t203 * t181;
t174 = t203 * t180 + t201 * t181;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t187, t188, 0, -t202 * t192 + t204 * t193, 0, 0, 0, 0, 0, 0, t204 * t178 - t202 * t185, t204 * t179 + t202 * t184, t204 * t186 - t202 * t189, t204 * t175 - t202 * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t176, t177, 0, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t188, t187, 0, -t204 * t192 - t202 * t193, 0, 0, 0, 0, 0, 0, -t202 * t178 - t204 * t185, -t202 * t179 + t204 * t184, -t202 * t186 - t204 * t189, -t202 * t175 - t204 * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, -qJDD(1), 0, t193, 0, 0, 0, 0, 0, 0, t178, t179, t186, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t206, 0, t192, 0, 0, 0, 0, 0, 0, t185, -t184, t189, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t176, t177, 0, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t191, t207, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t194, -t208, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t184, -t189, -t183;];
f_new_reg  = t1;
