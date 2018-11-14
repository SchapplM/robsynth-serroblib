% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:34
% EndTime: 2018-11-14 13:49:34
% DurationCPUTime: 0.02s
% Computational Cost: add. (24->12), mult. (44->26), div. (0->0), fcn. (8->2), ass. (0->12)
t43 = qJD(1) ^ 2;
t46 = t43 / 0.2e1;
t45 = qJ(2) * qJD(1);
t38 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t44 = t42 * t38 - t41 * t45;
t40 = -qJD(1) + qJD(3);
t39 = -qJD(1) * pkin(1) + qJD(2);
t36 = t41 * t38 + t42 * t45;
t35 = t40 * pkin(3) + t44;
t1 = [t46, 0, 0, -t39 * qJD(1), t43 * qJ(2), qJ(2) ^ 2 * t46 + t39 ^ 2 / 0.2e1, t40 ^ 2 / 0.2e1, t44 * t40, -t36 * t40, t36 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t1;
