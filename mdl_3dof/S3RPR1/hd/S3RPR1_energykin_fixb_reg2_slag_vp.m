% Calculate inertial parameters regressor of fixed base kinetic energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% T_reg [1x(3*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3RPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energykin_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_energykin_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:28
% EndTime: 2018-11-14 10:14:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (21->10), mult. (43->23), div. (0->0), fcn. (8->2), ass. (0->11)
t9 = qJD(1) ^ 2;
t6 = t9 / 0.2e1;
t10 = qJ(2) * qJD(1);
t8 = cos(qJ(3));
t7 = sin(qJ(3));
t5 = -qJD(1) + qJD(3);
t4 = -qJD(1) * pkin(1) + qJD(2);
t3 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t2 = t8 * t10 + t7 * t3;
t1 = -t7 * t10 + t8 * t3;
t11 = [0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, -t4 * qJD(1), 0, t9 * qJ(2), qJ(2) ^ 2 * t6 + t4 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t5 ^ 2 / 0.2e1, t1 * t5, -t2 * t5, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t11;
