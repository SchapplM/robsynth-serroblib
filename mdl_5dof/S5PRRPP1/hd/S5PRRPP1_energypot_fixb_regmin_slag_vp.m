% Calculate minimal parameter regressor of potential energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:38
% EndTime: 2021-01-15 15:22:38
% DurationCPUTime: 0.05s
% Computational Cost: add. (97->30), mult. (73->33), div. (0->0), fcn. (66->8), ass. (0->19)
t58 = cos(qJ(3));
t46 = t58 * pkin(3) + pkin(2);
t54 = pkin(7) + qJ(2);
t47 = sin(t54);
t49 = cos(t54);
t56 = -pkin(6) - qJ(4);
t63 = t47 * t46 + t49 * t56 + sin(pkin(7)) * pkin(1);
t57 = sin(qJ(3));
t62 = t57 * pkin(3) + pkin(5) + qJ(1);
t61 = t49 * t46 - t47 * t56 + cos(pkin(7)) * pkin(1);
t60 = g(1) * t49 + g(2) * t47;
t55 = qJ(3) + pkin(8);
t48 = sin(t55);
t50 = cos(t55);
t59 = pkin(4) * t50 + qJ(5) * t48;
t41 = g(1) * t47 - g(2) * t49;
t40 = -g(3) * t48 - t60 * t50;
t39 = -g(3) * t50 + t60 * t48;
t1 = [-g(3) * qJ(1), 0, -t60, t41, 0, 0, 0, 0, 0, -g(3) * t57 - t60 * t58, -g(3) * t58 + t60 * t57, t40, t39, -t41, -g(1) * t61 - g(2) * t63 - g(3) * t62, t40, -t41, -t39, -g(1) * (t59 * t49 + t61) - g(2) * (t59 * t47 + t63) - g(3) * (t48 * pkin(4) - t50 * qJ(5) + t62);];
U_reg = t1;
