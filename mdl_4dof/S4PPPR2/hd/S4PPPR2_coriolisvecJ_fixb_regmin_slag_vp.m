% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% taug_reg [4x6]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:55:35
% EndTime: 2018-11-14 13:55:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->5), mult. (36->9), div. (0->0), fcn. (26->4), ass. (0->6)
t6 = qJD(4) ^ 2;
t5 = cos(qJ(4));
t4 = sin(qJ(4));
t3 = cos(pkin(5));
t2 = sin(pkin(5));
t1 = [0, 0, 0, 0 (-t2 * t5 + t3 * t4) * t6 -(-t2 * t4 - t3 * t5) * t6; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t6 * t4, -t6 * t5; 0, 0, 0, 0, 0, 0;];
tauc_reg  = t1;
