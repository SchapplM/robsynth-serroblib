% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:29
% EndTime: 2018-11-14 14:09:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (24->12), mult. (56->27), div. (0->0), fcn. (26->4), ass. (0->13)
t45 = cos(qJ(2));
t39 = qJD(2) * pkin(2) + t45 * qJD(1);
t42 = sin(pkin(5));
t43 = cos(pkin(5));
t44 = sin(qJ(2));
t47 = qJD(1) * t44;
t37 = t42 * t39 + t43 * t47;
t46 = qJD(1) * qJD(2);
t36 = t43 * t39 - t42 * t47;
t41 = qJD(3) ^ 2 / 0.2e1;
t35 = qJD(2) * qJ(4) + t37;
t34 = -qJD(2) * pkin(3) + qJD(4) - t36;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t45 * t46, -t44 * t46, t37 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1 + t41, -t34 * qJD(2), t35 * qJD(2), t35 ^ 2 / 0.2e1 + t41 + t34 ^ 2 / 0.2e1;];
T_reg  = t1;
