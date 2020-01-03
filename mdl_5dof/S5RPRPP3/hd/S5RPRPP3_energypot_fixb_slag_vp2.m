% Calculate potential energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (125->50), mult. (135->39), div. (0->0), fcn. (104->6), ass. (0->20)
t54 = pkin(7) + qJ(3);
t51 = sin(t54);
t52 = cos(t54);
t79 = pkin(3) * t52 + qJ(4) * t51;
t78 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t77 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t76 = -m(5) - m(6);
t59 = cos(qJ(1));
t75 = t79 * t59;
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t73 = -m(3) * pkin(1) - t56 * mrSges(3,1) + t55 * mrSges(3,2) + t77 * t51 + t78 * t52 - mrSges(2,1);
t72 = -m(3) * qJ(2) - mrSges(5,1) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t70 = t55 * pkin(2) + pkin(5);
t49 = t56 * pkin(2) + pkin(1);
t45 = t59 * t49;
t57 = -pkin(6) - qJ(2);
t58 = sin(qJ(1));
t66 = -t58 * t57 + t45;
t1 = (-m(4) * t70 - t55 * mrSges(3,1) - t56 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t76 * (t51 * pkin(3) - qJ(4) * t52 + t70) + (-m(2) - m(3)) * pkin(5) - t77 * t52 + t78 * t51) * g(3) + (-mrSges(1,2) + (m(6) * pkin(4) - t72) * t59 + (-m(4) + t76) * (t58 * t49 + t59 * t57) + (t76 * t79 + t73) * t58) * g(2) + (-mrSges(1,1) - m(4) * t66 - m(5) * (t66 + t75) - m(6) * (t45 + t75) + t73 * t59 + (-m(6) * (pkin(4) - t57) + t72) * t58) * g(1);
U = t1;
