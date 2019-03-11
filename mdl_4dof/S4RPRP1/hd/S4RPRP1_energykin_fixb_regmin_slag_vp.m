% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:42
% EndTime: 2019-03-08 18:29:42
% DurationCPUTime: 0.02s
% Computational Cost: add. (33->13), mult. (74->29), div. (0->0), fcn. (26->4), ass. (0->14)
t52 = cos(pkin(6));
t47 = (pkin(1) * t52 + pkin(2)) * qJD(1);
t53 = sin(qJ(3));
t54 = cos(qJ(3));
t51 = sin(pkin(6));
t58 = pkin(1) * qJD(1) * t51;
t59 = t53 * t47 + t54 * t58;
t57 = t54 * t47 - t53 * t58;
t55 = qJD(1) ^ 2;
t50 = qJD(2) ^ 2 / 0.2e1;
t49 = qJD(1) + qJD(3);
t45 = t49 * qJ(4) + t59;
t44 = -t49 * pkin(3) + qJD(4) - t57;
t1 = [t55 / 0.2e1, 0, 0, t50 + (t51 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t55, t49 ^ 2 / 0.2e1, t57 * t49, -t59 * t49, -t44 * t49, t45 * t49, t45 ^ 2 / 0.2e1 + t50 + t44 ^ 2 / 0.2e1;];
T_reg  = t1;
