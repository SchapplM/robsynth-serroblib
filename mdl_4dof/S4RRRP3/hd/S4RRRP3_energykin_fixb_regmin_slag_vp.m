% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:14
% EndTime: 2019-12-31 17:14:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->17), mult. (110->47), div. (0->0), fcn. (40->4), ass. (0->18)
t63 = qJD(1) + qJD(2);
t62 = t63 ^ 2;
t75 = t62 / 0.2e1;
t64 = sin(qJ(3));
t74 = t63 * t64;
t66 = cos(qJ(3));
t73 = t63 * t66;
t72 = pkin(1) * qJD(1);
t71 = qJD(3) * t64;
t70 = qJD(3) * t66;
t69 = sin(qJ(2)) * t72;
t68 = cos(qJ(2)) * t72;
t61 = -t63 * pkin(2) - t68;
t60 = t63 * pkin(6) + t69;
t59 = qJD(3) * qJ(4) + t66 * t60;
t58 = -qJD(3) * pkin(3) + t64 * t60 + qJD(4);
t57 = -t68 + (-pkin(3) * t66 - qJ(4) * t64 - pkin(2)) * t63;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t75, t63 * t68, -t63 * t69, t64 ^ 2 * t75, t64 * t62 * t66, t63 * t71, t63 * t70, qJD(3) ^ 2 / 0.2e1, -t60 * t71 - t61 * t73, -t60 * t70 + t61 * t74, -t58 * qJD(3) - t57 * t73, (t58 * t64 + t59 * t66) * t63, t59 * qJD(3) - t57 * t74, t59 ^ 2 / 0.2e1 + t57 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1;];
T_reg = t1;
