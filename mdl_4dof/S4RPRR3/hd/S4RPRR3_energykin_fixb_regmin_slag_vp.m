% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (44->22), mult. (134->60), div. (0->0), fcn. (69->6), ass. (0->23)
t81 = qJD(1) ^ 2;
t88 = t81 / 0.2e1;
t75 = sin(pkin(7));
t69 = (pkin(1) * t75 + pkin(5)) * qJD(1);
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t87 = t78 * qJD(2) + t80 * t69;
t86 = qJD(1) * t78;
t85 = qJD(1) * t80;
t84 = qJD(1) * qJD(3);
t76 = cos(pkin(7));
t83 = -pkin(1) * t76 - pkin(2);
t79 = cos(qJ(4));
t77 = sin(qJ(4));
t74 = qJD(3) + qJD(4);
t73 = t80 * qJD(2);
t70 = t83 * qJD(1);
t67 = (-pkin(3) * t80 + t83) * qJD(1);
t66 = (t77 * t80 + t78 * t79) * qJD(1);
t65 = t77 * t86 - t79 * t85;
t64 = pkin(6) * t85 + t87;
t63 = qJD(3) * pkin(3) + t73 + (-pkin(6) * qJD(1) - t69) * t78;
t1 = [t88, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t75 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t81, t78 ^ 2 * t88, t78 * t81 * t80, t78 * t84, t80 * t84, qJD(3) ^ 2 / 0.2e1, -t70 * t85 + (-t78 * t69 + t73) * qJD(3), -t87 * qJD(3) + t70 * t86, t66 ^ 2 / 0.2e1, -t66 * t65, t66 * t74, -t65 * t74, t74 ^ 2 / 0.2e1, t67 * t65 + (t79 * t63 - t77 * t64) * t74, t67 * t66 - (t77 * t63 + t79 * t64) * t74;];
T_reg = t1;
