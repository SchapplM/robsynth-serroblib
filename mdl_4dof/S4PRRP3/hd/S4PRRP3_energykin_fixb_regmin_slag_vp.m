% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:29
% EndTime: 2018-11-14 14:12:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->9), mult. (39->25), div. (0->0), fcn. (18->4), ass. (0->12)
t36 = sin(qJ(2));
t41 = qJD(1) * t36;
t40 = qJD(1) * qJD(2);
t38 = cos(qJ(2));
t33 = qJD(2) * pkin(2) + t38 * qJD(1);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t39 = t37 * t33 - t35 * t41;
t34 = qJD(2) + qJD(3);
t31 = t35 * t33 + t37 * t41;
t30 = t34 * pkin(3) + t39;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t38 * t40, -t36 * t40, t34 ^ 2 / 0.2e1, t39 * t34, -t31 * t34, t31 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t1;
