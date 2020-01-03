% Calculate potential energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:34
% EndTime: 2019-12-31 17:39:34
% DurationCPUTime: 0.14s
% Computational Cost: add. (135->56), mult. (108->62), div. (0->0), fcn. (98->8), ass. (0->22)
t61 = rSges(6,3) + pkin(7);
t60 = cos(qJ(4));
t59 = sin(qJ(4));
t58 = pkin(5) + qJ(1);
t46 = sin(pkin(8));
t44 = t46 * pkin(1);
t56 = pkin(8) + qJ(2);
t52 = sin(t56);
t57 = t52 * pkin(2) + t44;
t55 = -pkin(6) + t58;
t43 = cos(t56);
t47 = cos(pkin(8));
t45 = t47 * pkin(1);
t54 = t43 * pkin(2) + t52 * qJ(3) + t45;
t53 = t43 * pkin(3) + t54;
t51 = t52 * pkin(3) - t43 * qJ(3) + t57;
t48 = sin(qJ(5));
t49 = cos(qJ(5));
t50 = -rSges(6,1) * t49 + rSges(6,2) * t48 - pkin(4);
t34 = t43 * t59 - t52 * t60;
t33 = -t43 * t60 - t52 * t59;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t47 * rSges(2,1) - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) + t47 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t43 * rSges(3,1) - rSges(3,2) * t52 + t45) + g(2) * (rSges(3,1) * t52 + t43 * rSges(3,2) + t44) + g(3) * (rSges(3,3) + t58)) - m(4) * (g(1) * (t43 * rSges(4,1) + rSges(4,3) * t52 + t54) + g(2) * (t52 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t43 + t57) + g(3) * (rSges(4,2) + t58)) - m(5) * (g(1) * (-t33 * rSges(5,1) - t34 * rSges(5,2) + t53) + g(2) * (-t34 * rSges(5,1) + t33 * rSges(5,2) + t51) + g(3) * (-rSges(5,3) + t55)) - m(6) * (g(1) * t53 + g(2) * t51 + g(3) * (-t48 * rSges(6,1) - t49 * rSges(6,2) + t55) + (g(1) * t61 + g(2) * t50) * t34 + (g(1) * t50 - g(2) * t61) * t33);
U = t1;
