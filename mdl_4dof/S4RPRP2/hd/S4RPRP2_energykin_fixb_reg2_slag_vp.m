% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:34
% EndTime: 2018-11-14 13:49:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->14), mult. (75->27), div. (0->0), fcn. (16->2), ass. (0->15)
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t14 = qJ(2) * qJD(1);
t6 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t4 = t11 * t6 + t12 * t14;
t9 = -qJD(1) + qJD(3);
t15 = t4 * t9;
t13 = qJD(1) ^ 2;
t10 = t13 / 0.2e1;
t3 = -t11 * t14 + t12 * t6;
t8 = -qJD(1) * pkin(1) + qJD(2);
t7 = t9 ^ 2 / 0.2e1;
t2 = t4 ^ 2 / 0.2e1;
t1 = t9 * pkin(3) + t3;
t5 = [0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, -t8 * qJD(1), 0, t13 * qJ(2), qJ(2) ^ 2 * t10 + t8 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t7, t3 * t9, -t15, 0, t2 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t7, t1 * t9, -t15, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t5;
