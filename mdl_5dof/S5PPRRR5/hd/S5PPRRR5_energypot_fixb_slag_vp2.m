% Calculate potential energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:36
% EndTime: 2019-12-31 17:35:36
% DurationCPUTime: 0.24s
% Computational Cost: add. (120->56), mult. (133->49), div. (0->0), fcn. (120->8), ass. (0->27)
t77 = -m(3) - m(4);
t76 = -m(5) - m(6);
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t75 = m(6) * pkin(4) + t60 * mrSges(6,1) - t58 * mrSges(6,2) + mrSges(5,1);
t74 = mrSges(2,2) - mrSges(3,3);
t73 = -m(4) * pkin(2) - mrSges(2,1) - mrSges(3,1);
t72 = m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3);
t56 = sin(pkin(8));
t59 = sin(qJ(3));
t71 = t56 * t59;
t57 = cos(pkin(8));
t70 = t57 * t59;
t69 = -pkin(5) + qJ(1);
t68 = t57 * pkin(1) + t56 * qJ(2);
t67 = qJ(3) + qJ(4);
t66 = cos(t67);
t65 = sin(t67);
t61 = cos(qJ(3));
t53 = t56 * pkin(1);
t51 = t61 * pkin(3) + pkin(2);
t45 = t56 * t51;
t44 = t56 * t61 - t70;
t43 = -t57 * t61 - t71;
t42 = -t56 * t66 + t57 * t65;
t41 = -t56 * t65 - t57 * t66;
t1 = (-m(4) * t69 + t58 * mrSges(6,1) + t60 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + mrSges(5,3) + t76 * (-pkin(6) + t69) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-mrSges(1,2) - t44 * mrSges(4,1) - t43 * mrSges(4,2) - m(5) * (t45 + t53) - m(6) * (-pkin(3) * t70 + t45) + (-m(6) + t77) * (-t57 * qJ(2) + t53) + (-m(5) * (-pkin(3) * t59 - qJ(2)) - t74) * t57 + t73 * t56 + t75 * t42 + t72 * t41) * g(2) + (t43 * mrSges(4,1) - t44 * mrSges(4,2) - mrSges(1,1) + t77 * t68 + t76 * (pkin(3) * t71 + t57 * t51 + t68) + t73 * t57 + t74 * t56 - t72 * t42 + t75 * t41) * g(1);
U = t1;
