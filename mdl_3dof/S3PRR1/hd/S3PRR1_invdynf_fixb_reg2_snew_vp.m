% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S3PRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
%
% Output:
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3PRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:25:02
% EndTime: 2019-05-04 18:25:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (290->43), mult. (358->36), div. (0->0), fcn. (260->4), ass. (0->28)
t203 = qJD(2) + qJD(3);
t201 = t203 ^ 2;
t202 = qJDD(2) + qJDD(3);
t205 = sin(qJ(3));
t207 = cos(qJ(3));
t189 = t205 * t201 - t207 * t202;
t206 = sin(qJ(2));
t208 = cos(qJ(2));
t210 = -t207 * t201 - t205 * t202;
t213 = t206 * t189 + t208 * t210;
t181 = t208 * t189 - t206 * t210;
t204 = -g(2) + qJDD(1);
t193 = t206 * g(1) + t208 * t204;
t194 = -t208 * g(1) + t206 * t204;
t209 = qJD(2) ^ 2;
t196 = t208 * qJDD(2) - t206 * t209;
t195 = -t206 * qJDD(2) - t208 * t209;
t192 = -t209 * pkin(2) + t194;
t191 = qJDD(2) * pkin(2) + t193;
t186 = -t206 * t193 + t208 * t194;
t185 = t208 * t193 + t206 * t194;
t184 = t205 * t191 + t207 * t192;
t183 = t207 * t191 - t205 * t192;
t178 = -t205 * t183 + t207 * t184;
t177 = t207 * t183 + t205 * t184;
t176 = -t206 * t177 + t208 * t178;
t175 = t208 * t177 + t206 * t178;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t195, -t196, 0, t186, 0, 0, 0, 0, 0, 0, t213, t181, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, t196, t195, 0, t185, 0, 0, 0, 0, 0, 0, -t181, t213, 0, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t195, -t196, 0, t186, 0, 0, 0, 0, 0, 0, t213, t181, 0, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, t196, t195, 0, t185, 0, 0, 0, 0, 0, 0, -t181, t213, 0, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, -qJDD(2), 0, t194, 0, 0, 0, 0, 0, 0, t210, t189, 0, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t209, 0, t193, 0, 0, 0, 0, 0, 0, -t189, t210, 0, t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, -t202, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t201, 0, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
f_new_reg  = t1;
