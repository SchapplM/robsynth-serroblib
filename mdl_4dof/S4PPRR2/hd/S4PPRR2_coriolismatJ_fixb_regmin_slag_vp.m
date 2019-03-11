% Calculate minimal parameter regressor of coriolis matrix for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x8]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (70->10), mult. (170->21), div. (0->0), fcn. (196->6), ass. (0->15)
t27 = qJD(3) + qJD(4);
t12 = sin(pkin(6));
t14 = sin(qJ(3));
t17 = cos(pkin(6));
t24 = cos(qJ(3));
t10 = t24 * t12 + t14 * t17;
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t9 = t14 * t12 - t24 * t17;
t26 = t27 * (-t15 * t10 + t13 * t9);
t25 = t27 * (t13 * t10 + t15 * t9);
t23 = pkin(3) * qJD(4);
t22 = qJD(3) * pkin(3);
t16 = pkin(3) * t27;
t1 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t10 * qJD(3), t9 * qJD(3), 0, t26, t25; 0, 0, 0, 0, 0, 0, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t13 * t23, -t15 * t23; 0, 0, 0, 0, 0, 0, -t13 * t16, -t15 * t16; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t13 * t22, t15 * t22; 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
