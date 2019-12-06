% Calculate potential energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:29
% EndTime: 2019-12-05 15:37:30
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->81), mult. (195->101), div. (0->0), fcn. (187->8), ass. (0->31)
t88 = rSges(6,1) + pkin(4);
t64 = sin(pkin(7));
t66 = cos(pkin(7));
t87 = g(1) * t66 + g(2) * t64;
t86 = rSges(4,3) + qJ(3);
t85 = rSges(6,3) + qJ(5);
t68 = sin(qJ(2));
t82 = rSges(3,2) * t68;
t63 = sin(pkin(8));
t81 = t64 * t63;
t69 = cos(qJ(2));
t80 = t64 * t69;
t79 = t66 * t63;
t78 = t66 * t69;
t75 = t66 * pkin(1) + t64 * pkin(5);
t65 = cos(pkin(8));
t55 = pkin(3) * t65 + pkin(2);
t67 = -pkin(6) - qJ(3);
t73 = t68 * t55 + t69 * t67 + qJ(1);
t60 = t64 * pkin(1);
t72 = -t66 * pkin(5) + t60;
t71 = pkin(3) * t81 + t55 * t78 + t75;
t70 = -pkin(3) * t79 + t55 * t80 + t72;
t62 = pkin(8) + qJ(4);
t58 = cos(t62);
t57 = sin(t62);
t49 = t64 * t57 + t58 * t78;
t48 = t57 * t78 - t64 * t58;
t47 = -t66 * t57 + t58 * t80;
t46 = t57 * t80 + t66 * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t66 - rSges(2,2) * t64) + g(2) * (rSges(2,1) * t64 + rSges(2,2) * t66) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t64 * rSges(3,3) + t75) + g(2) * (rSges(3,1) * t80 - t64 * t82 + t60) + g(3) * (t68 * rSges(3,1) + rSges(3,2) * t69 + qJ(1)) + (g(1) * (rSges(3,1) * t69 - t82) + g(2) * (-rSges(3,3) - pkin(5))) * t66) - m(4) * (g(1) * (pkin(2) * t78 + (t65 * t78 + t81) * rSges(4,1) + (-t63 * t78 + t64 * t65) * rSges(4,2) + t75) + g(2) * (pkin(2) * t80 + (t65 * t80 - t79) * rSges(4,1) + (-t63 * t80 - t66 * t65) * rSges(4,2) + t72) + g(3) * (-t86 * t69 + qJ(1)) + (g(3) * (rSges(4,1) * t65 - rSges(4,2) * t63 + pkin(2)) + t87 * t86) * t68) - m(5) * (g(1) * (rSges(5,1) * t49 - rSges(5,2) * t48 + t71) + g(2) * (t47 * rSges(5,1) - rSges(5,2) * t46 + t70) + g(3) * (-rSges(5,3) * t69 + t73) + (g(3) * (rSges(5,1) * t58 - rSges(5,2) * t57) + t87 * (rSges(5,3) - t67)) * t68) - m(6) * (g(1) * (t85 * t48 + t88 * t49 + t71) + g(2) * (t85 * t46 + t88 * t47 + t70) + g(3) * (-rSges(6,2) * t69 + t73) + (g(3) * (t85 * t57 + t88 * t58) + t87 * (rSges(6,2) - t67)) * t68);
U = t1;
