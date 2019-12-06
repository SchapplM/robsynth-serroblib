% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:12
% EndTime: 2019-12-05 15:45:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (65->17), mult. (133->49), div. (0->0), fcn. (81->8), ass. (0->23)
t84 = qJD(2) + qJD(4);
t83 = t84 ^ 2;
t99 = t83 / 0.2e1;
t92 = cos(qJ(2));
t82 = qJD(2) * pkin(2) + t92 * qJD(1);
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t89 = sin(qJ(2));
t96 = qJD(1) * t89;
t79 = t86 * t82 - t85 * t96;
t78 = qJD(2) * pkin(3) + t79;
t80 = t85 * t82 + t86 * t96;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t93 = t91 * t78 - t88 * t80;
t98 = (-t84 * pkin(4) - t93) * t84;
t97 = t88 * t78 + t91 * t80;
t95 = qJD(5) * t84;
t94 = qJD(1) * qJD(2);
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t75 = t84 * pkin(7) + t97;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t92 * t94, -t89 * t94, t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t99, t93 * t84, -t97 * t84, t87 ^ 2 * t99, t87 * t83 * t90, t87 * t95, t90 * t95, qJD(5) ^ 2 / 0.2e1, (t90 * qJD(3) - t87 * t75) * qJD(5) - t90 * t98, -(t87 * qJD(3) + t90 * t75) * qJD(5) + t87 * t98;];
T_reg = t1;
