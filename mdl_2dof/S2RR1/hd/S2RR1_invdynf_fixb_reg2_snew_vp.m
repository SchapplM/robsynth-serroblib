% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S2RR1
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
% Datum: 2019-05-04 18:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S2RR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynf_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynf_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynf_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:17:53
% EndTime: 2019-05-04 18:17:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (127->41), mult. (306->60), div. (0->0), fcn. (194->4), ass. (0->36)
t205 = sin(qJ(1));
t207 = cos(qJ(1));
t194 = t205 * g(1) + t207 * g(3);
t204 = sin(qJ(2));
t202 = t204 ^ 2;
t206 = cos(qJ(2));
t203 = t206 ^ 2;
t213 = -t202 - t203;
t212 = t204 * qJDD(1);
t211 = t206 * qJDD(1);
t210 = -2 * qJD(1) * qJD(2);
t195 = t207 * g(1) - t205 * g(3);
t209 = qJD(1) ^ 2;
t208 = qJD(2) ^ 2;
t198 = t206 * t209 * t204;
t197 = -t203 * t209 - t208;
t196 = -t202 * t209 - t208;
t193 = -qJDD(2) + t198;
t192 = qJDD(2) + t198;
t191 = t213 * t209;
t190 = -t207 * qJDD(1) + t205 * t209;
t189 = t205 * qJDD(1) + t207 * t209;
t188 = t213 * qJDD(1);
t187 = t204 * t210 + t211;
t186 = t206 * t210 - t212;
t185 = -(t209 * pkin(1)) + t195;
t184 = -qJDD(1) * pkin(1) + t194;
t183 = t204 * g(2) + t206 * t184;
t182 = t206 * g(2) - t204 * t184;
t181 = t206 * t193 - t204 * t196;
t180 = -t204 * t192 + t206 * t197;
t179 = -t204 * t193 - t206 * t196;
t178 = -t206 * t192 - t204 * t197;
t177 = -t204 * t182 + t206 * t183;
t176 = -t206 * t182 - t204 * t183;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t190, t189, 0, -t205 * t194 - t207 * t195, 0, 0, 0, 0, 0, 0, -t205 * t180 - t207 * t187, -t205 * t181 - t207 * t186, -t205 * t188 - t207 * t191, -t205 * t177 - t207 * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t178, t179, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t189, -t190, 0, -t207 * t194 + t205 * t195, 0, 0, 0, 0, 0, 0, -t207 * t180 + t205 * t187, -t207 * t181 + t205 * t186, -t207 * t188 + t205 * t191, -t207 * t177 + t205 * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, -qJDD(1), 0, t194, 0, 0, 0, 0, 0, 0, t180, t181, t188, t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t209, 0, t195, 0, 0, 0, 0, 0, 0, t187, t186, t191, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t178, t179, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t193, -t211, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t196, t212, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t186, t191, t185;];
f_new_reg  = t1;
