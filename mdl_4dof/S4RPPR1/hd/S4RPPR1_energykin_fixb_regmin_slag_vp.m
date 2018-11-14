% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:45
% EndTime: 2018-11-14 13:46:45
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->14), mult. (58->29), div. (0->0), fcn. (14->4), ass. (0->12)
t47 = cos(pkin(6));
t52 = -pkin(1) * t47 - pkin(2);
t50 = qJD(1) ^ 2;
t49 = cos(qJ(4));
t48 = sin(qJ(4));
t46 = sin(pkin(6));
t45 = qJD(2) ^ 2 / 0.2e1;
t44 = -qJD(1) + qJD(4);
t43 = (pkin(1) * t46 + qJ(3)) * qJD(1);
t42 = t52 * qJD(1) + qJD(3);
t41 = qJD(3) + (-pkin(3) + t52) * qJD(1);
t1 = [t50 / 0.2e1, 0, 0, t45 + (t46 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t50, -t42 * qJD(1), t43 * qJD(1), t43 ^ 2 / 0.2e1 + t45 + t42 ^ 2 / 0.2e1, t44 ^ 2 / 0.2e1 (t49 * t41 - t48 * t43) * t44 -(t48 * t41 + t49 * t43) * t44;];
T_reg  = t1;
