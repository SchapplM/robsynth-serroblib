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
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (127->41), mult. (306->60), div. (0->0), fcn. (194->4), ass. (0->36)
t202 = sin(qJ(2));
t200 = t202 ^ 2;
t204 = cos(qJ(2));
t201 = t204 ^ 2;
t211 = -t200 - t201;
t210 = t202 * qJDD(1);
t209 = t204 * qJDD(1);
t208 = -2 * qJD(1) * qJD(2);
t203 = sin(qJ(1));
t205 = cos(qJ(1));
t195 = -t205 * g(1) + t203 * g(3);
t194 = -t203 * g(1) - t205 * g(3);
t207 = qJD(1) ^ 2;
t206 = qJD(2) ^ 2;
t198 = t204 * t207 * t202;
t197 = -t201 * t207 - t206;
t196 = -t200 * t207 - t206;
t193 = -qJDD(2) + t198;
t192 = qJDD(2) + t198;
t191 = t211 * t207;
t190 = t205 * qJDD(1) - t203 * t207;
t189 = -t203 * qJDD(1) - t205 * t207;
t188 = t211 * qJDD(1);
t187 = t202 * t208 + t209;
t186 = t204 * t208 - t210;
t185 = -(t207 * pkin(1)) + t195;
t184 = -qJDD(1) * pkin(1) + t194;
t183 = t202 * g(2) + t204 * t184;
t182 = t204 * g(2) - t202 * t184;
t181 = t204 * t193 - t202 * t196;
t180 = -t202 * t192 + t204 * t197;
t179 = -t202 * t193 - t204 * t196;
t178 = -t204 * t192 - t202 * t197;
t177 = -t202 * t182 + t204 * t183;
t176 = -t204 * t182 - t202 * t183;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t190, t189, 0, t203 * t194 + t205 * t195, 0, 0, 0, 0, 0, 0, t203 * t180 + t205 * t187, t203 * t181 + t205 * t186, t203 * t188 + t205 * t191, t203 * t177 + t205 * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t178, t179, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t189, -t190, 0, t205 * t194 - t203 * t195, 0, 0, 0, 0, 0, 0, t205 * t180 - t203 * t187, t205 * t181 - t203 * t186, t205 * t188 - t203 * t191, t205 * t177 - t203 * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -qJDD(1), 0, t194, 0, 0, 0, 0, 0, 0, t180, t181, t188, t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t207, 0, t195, 0, 0, 0, 0, 0, 0, t187, t186, t191, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t178, t179, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t193, -t209, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t196, t210, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t186, t191, t185;];
f_new_reg = t1;
