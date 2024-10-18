% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:55
% EndTime: 2024-09-27 21:45:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (152->35), mult. (460->78), div. (0->0), fcn. (356->8), ass. (0->35)
t44 = sin(pkin(5));
t45 = cos(pkin(5));
t57 = t45 * qJD(2);
t64 = pkin(2) * t57 + qJD(1) * t44;
t63 = cos(qJ(4));
t51 = qJD(2) ^ 2;
t62 = t44 ^ 2 * t51;
t41 = qJD(3) + t57;
t48 = sin(qJ(3));
t58 = qJD(2) * t44;
t54 = t48 * t58;
t50 = cos(qJ(3));
t60 = t64 * t50;
t25 = t41 * pkin(3) + (-pkin(7) - pkin(8)) * t54 + t60;
t53 = t50 * t58;
t56 = pkin(7) * t53 + t64 * t48;
t27 = pkin(8) * t53 + t56;
t47 = sin(qJ(4));
t61 = t47 * t25 + t63 * t27;
t52 = t63 * t25 - t47 * t27;
t35 = qJD(4) + t41;
t42 = t45 * qJD(1);
t31 = t42 + (-pkin(3) * t50 - pkin(2)) * t58;
t49 = cos(qJ(5));
t46 = sin(qJ(5));
t34 = qJD(5) + t35;
t32 = -pkin(2) * t58 + t42;
t30 = (t47 * t50 + t63 * t48) * t58;
t29 = t47 * t54 - t63 * t53;
t24 = t29 * pkin(4) + t31;
t21 = -t46 * t29 + t49 * t30;
t20 = t49 * t29 + t46 * t30;
t19 = -t29 * pkin(9) + t61;
t18 = t35 * pkin(4) - t30 * pkin(9) + t52;
t1 = [qJD(1) ^ 2 / 0.2e1, t51 / 0.2e1, 0, 0, t48 ^ 2 * t62 / 0.2e1, t48 * t50 * t62, t41 * t54, t41 * t53, t41 ^ 2 / 0.2e1, (-pkin(7) * t54 + t60) * t41 - t32 * t53, t32 * t54 - t56 * t41, t30 ^ 2 / 0.2e1, -t30 * t29, t30 * t35, -t29 * t35, t35 ^ 2 / 0.2e1, t31 * t29 + t52 * t35, t31 * t30 - t61 * t35, t21 ^ 2 / 0.2e1, -t21 * t20, t21 * t34, -t20 * t34, t34 ^ 2 / 0.2e1, (t49 * t18 - t46 * t19) * t34 + t24 * t20, -(t46 * t18 + t49 * t19) * t34 + t24 * t21;];
T_reg = t1;
