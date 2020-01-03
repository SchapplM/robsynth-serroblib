% Calculate inertial parameters regressor of coriolis matrix for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% cmat_reg [(2*2)x(2*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S2RR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.05s
% Computational Cost: add. (8->7), mult. (27->10), div. (0->0), fcn. (22->2), ass. (0->9)
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t1 = -t3 ^ 2 + t4 ^ 2;
t8 = t1 * qJD(1);
t7 = t3 * qJD(2);
t6 = t4 * qJD(2);
t5 = t3 * t6;
t2 = t3 * t4 * qJD(1);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t1 * qJD(2), 0, -t5, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t8, -t6, -t2, t7, 0, pkin(1) * t6, -pkin(1) * t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t8, 0, t2, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
