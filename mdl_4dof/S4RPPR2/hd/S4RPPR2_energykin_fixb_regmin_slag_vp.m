% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% T_reg [1x12]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:31
% EndTime: 2018-11-14 13:47:31
% DurationCPUTime: 0.03s
% Computational Cost: add. (39->16), mult. (71->33), div. (0->0), fcn. (20->4), ass. (0->15)
t51 = qJD(1) ^ 2;
t53 = t51 / 0.2e1;
t52 = qJ(2) * qJD(1);
t50 = cos(qJ(4));
t49 = sin(qJ(4));
t48 = cos(pkin(6));
t47 = sin(pkin(6));
t46 = -qJD(1) + qJD(4);
t45 = -qJD(1) * pkin(1) + qJD(2);
t44 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t43 = t48 * t44;
t42 = t47 * t44 + t48 * t52;
t41 = -t47 * t52 + t43;
t40 = t43 + (-qJ(2) * t47 - pkin(3)) * qJD(1);
t1 = [t53, 0, 0, -t45 * qJD(1), t51 * qJ(2), qJ(2) ^ 2 * t53 + t45 ^ 2 / 0.2e1, -t41 * qJD(1), t42 * qJD(1), t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t46 ^ 2 / 0.2e1 (t50 * t40 - t49 * t42) * t46 -(t49 * t40 + t50 * t42) * t46;];
T_reg  = t1;
