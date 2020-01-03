% Calculate potential energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:48
% EndTime: 2019-12-31 20:54:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (139->48), mult. (124->38), div. (0->0), fcn. (93->8), ass. (0->21)
t71 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t70 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t69 = -m(5) - m(6);
t51 = qJ(2) + qJ(3);
t45 = pkin(8) + t51;
t41 = sin(t45);
t42 = cos(t45);
t54 = cos(qJ(2));
t44 = t54 * pkin(2) + pkin(1);
t46 = sin(t51);
t47 = cos(t51);
t52 = sin(qJ(2));
t68 = -m(3) * pkin(1) - m(4) * t44 - t54 * mrSges(3,1) - t47 * mrSges(4,1) + t52 * mrSges(3,2) + t46 * mrSges(4,2) - t70 * t41 + t71 * t42 - mrSges(2,1);
t56 = -pkin(7) - pkin(6);
t67 = -m(3) * pkin(6) + m(4) * t56 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t66 = t52 * pkin(2) + pkin(5);
t55 = cos(qJ(1));
t53 = sin(qJ(1));
t50 = -qJ(4) + t56;
t39 = pkin(3) * t47 + t44;
t1 = (-m(4) * t66 - t52 * mrSges(3,1) - t46 * mrSges(4,1) - t54 * mrSges(3,2) - t47 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t69 * (pkin(3) * t46 + t66) + t70 * t42 + t71 * t41 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + t69 * (t53 * t39 + t55 * t50) - t67 * t55 + t68 * t53) * g(2) + (-mrSges(1,1) + t69 * (t55 * t39 - t53 * t50) + t68 * t55 + t67 * t53) * g(1);
U = t1;
