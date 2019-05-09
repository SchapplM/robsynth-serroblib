% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPPR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:40:49
% EndTime: 2019-05-04 18:40:50
% DurationCPUTime: 0.31s
% Computational Cost: add. (345->41), mult. (370->29), div. (0->0), fcn. (340->4), ass. (0->23)
t209 = sin(qJ(4));
t210 = cos(qJ(4));
t211 = qJD(4) ^ 2;
t201 = t210 * qJDD(4) - t209 * t211;
t202 = -t209 * qJDD(4) - t210 * t211;
t207 = sin(pkin(5));
t208 = cos(pkin(5));
t196 = t208 * t201 + t207 * t202;
t212 = -t207 * t201 + t208 * t202;
t206 = g(1) - qJDD(2);
t205 = -g(2) + qJDD(1);
t204 = g(3) + qJDD(3);
t200 = t208 * t205 - t207 * t206;
t199 = -t207 * t205 - t208 * t206;
t194 = t209 * t199 + t210 * t200;
t193 = t210 * t199 - t209 * t200;
t192 = -t207 * t199 + t208 * t200;
t191 = t208 * t199 + t207 * t200;
t190 = -t209 * t193 + t210 * t194;
t189 = t210 * t193 + t209 * t194;
t188 = -t207 * t189 + t208 * t190;
t187 = t208 * t189 + t207 * t190;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, t196, t212, 0, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, t212, -t196, 0, t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, 0, 0, 0, 0, 0, 0, -t196, -t212, 0, -t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, t212, -t196, 0, t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, t212, -t196, 0, t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, t196, t212, 0, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, 0, 0, 0, 0, 0, 0, t202, -t201, 0, t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, 0, 0, 0, 0, 0, 0, t201, t202, 0, t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, -qJDD(4), 0, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t211, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204;];
f_new_reg  = t1;
