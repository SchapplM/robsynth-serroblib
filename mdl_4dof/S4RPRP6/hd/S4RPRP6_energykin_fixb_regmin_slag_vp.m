% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP6
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
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:39
% EndTime: 2021-01-15 10:27:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (43->17), mult. (98->44), div. (0->0), fcn. (28->2), ass. (0->15)
t10 = qJD(1) ^ 2;
t16 = t10 / 0.2e1;
t8 = sin(qJ(3));
t5 = qJD(4) + (pkin(3) * t8 + qJ(2)) * qJD(1);
t15 = qJD(1) * t5;
t6 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t14 = qJD(3) * t6;
t13 = t10 * qJ(2);
t12 = qJD(1) * qJD(3);
t11 = -qJ(4) * qJD(1) + t6;
t9 = cos(qJ(3));
t7 = -qJD(1) * pkin(1) + qJD(2);
t4 = t11 * t8;
t3 = qJD(3) * pkin(3) + t11 * t9;
t1 = [t16, 0, 0, t7 * qJD(1), t13, qJ(2) ^ 2 * t16 + t7 ^ 2 / 0.2e1, t9 ^ 2 * t16, -t9 * t10 * t8, t9 * t12, -t8 * t12, qJD(3) ^ 2 / 0.2e1, t8 * t13 + t9 * t14, t9 * t13 - t8 * t14, t3 * qJD(3) + t8 * t15, -t4 * qJD(3) + t9 * t15, (-t3 * t9 - t4 * t8) * qJD(1), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t1;
