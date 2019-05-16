% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S3PPP1
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,theta1]';
%
% Output:
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3PPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:20:10
% EndTime: 2019-05-04 18:20:10
% DurationCPUTime: 0.13s
% Computational Cost: add. (56->21), mult. (60->12), div. (0->0), fcn. (60->2), ass. (0->12)
t79 = sin(pkin(3));
t80 = cos(pkin(3));
t75 = t79 * g(1) - t80 * g(2);
t73 = -qJDD(2) + t75;
t81 = t79 * t73;
t76 = t80 * g(1) + t79 * g(2);
t78 = -g(3) + qJDD(1);
t74 = -qJDD(3) + t76;
t72 = t80 * t76;
t71 = t79 * t76;
t70 = t80 * t73;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t75 - t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 - t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t74 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t75 - t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 + t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t74 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74;];
f_new_reg  = t1;
