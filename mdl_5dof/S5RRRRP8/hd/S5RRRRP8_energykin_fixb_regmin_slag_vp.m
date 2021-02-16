% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:19
% EndTime: 2021-01-16 00:21:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (263->38), mult. (586->83), div. (0->0), fcn. (389->6), ass. (0->34)
t54 = qJD(1) ^ 2;
t66 = t54 / 0.2e1;
t65 = cos(qJ(4));
t53 = cos(qJ(2));
t64 = t53 * t54;
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t51 = sin(qJ(2));
t61 = qJD(1) * t51;
t39 = t50 * qJD(2) + t52 * t61;
t60 = t53 * qJD(1);
t45 = -qJD(3) + t60;
t37 = (-pkin(2) * t53 - pkin(7) * t51 - pkin(1)) * qJD(1);
t42 = pkin(6) * t60 + qJD(2) * pkin(7);
t55 = t52 * t37 - t50 * t42;
t28 = -t45 * pkin(3) - t39 * pkin(8) + t55;
t38 = -t52 * qJD(2) + t50 * t61;
t62 = t50 * t37 + t52 * t42;
t30 = -t38 * pkin(8) + t62;
t49 = sin(qJ(4));
t63 = t49 * t28 + t65 * t30;
t59 = qJD(1) * qJD(2);
t58 = t51 * t59;
t57 = t53 * t59;
t56 = t65 * t28 - t49 * t30;
t41 = -qJD(2) * pkin(2) + pkin(6) * t61;
t33 = t38 * pkin(3) + t41;
t43 = -qJD(4) + t45;
t32 = -t49 * t38 + t65 * t39;
t31 = t65 * t38 + t49 * t39;
t25 = t31 * pkin(4) + qJD(5) + t33;
t24 = -t31 * qJ(5) + t63;
t23 = -t43 * pkin(4) - t32 * qJ(5) + t56;
t1 = [t66, 0, 0, t51 ^ 2 * t66, t51 * t64, t58, t57, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(6) * t58, -t54 * pkin(1) * t51 - pkin(6) * t57, t39 ^ 2 / 0.2e1, -t39 * t38, -t39 * t45, t38 * t45, t45 ^ 2 / 0.2e1, t41 * t38 - t55 * t45, t41 * t39 + t62 * t45, t32 ^ 2 / 0.2e1, -t32 * t31, -t32 * t43, t31 * t43, t43 ^ 2 / 0.2e1, t33 * t31 - t56 * t43, t33 * t32 + t63 * t43, -t23 * t43 + t25 * t31, t24 * t43 + t25 * t32, -t23 * t32 - t24 * t31, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1;];
T_reg = t1;
