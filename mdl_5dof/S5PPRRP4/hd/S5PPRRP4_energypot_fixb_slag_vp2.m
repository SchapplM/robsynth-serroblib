% Calculate potential energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:23
% EndTime: 2019-12-31 17:34:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (100->44), mult. (148->35), div. (0->0), fcn. (141->6), ass. (0->19)
t68 = m(5) + m(6);
t67 = mrSges(5,2) + mrSges(6,2);
t66 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t65 = -mrSges(2,1) - mrSges(3,1);
t64 = mrSges(2,2) - mrSges(3,3);
t63 = -m(4) - t68;
t50 = sin(qJ(4));
t51 = cos(qJ(4));
t62 = t68 * pkin(3) - t67 * t50 + t66 * t51 + mrSges(4,1);
t61 = m(5) * pkin(6) - m(6) * (-qJ(5) - pkin(6)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t60 = cos(qJ(3));
t59 = sin(qJ(3));
t47 = cos(pkin(7));
t57 = sin(pkin(7));
t58 = t47 * pkin(1) + t57 * qJ(2);
t55 = t57 * pkin(1) - t47 * qJ(2);
t37 = t47 * t59 - t57 * t60;
t36 = -t47 * t60 - t57 * t59;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t67 * t51 + t66 * t50 + t63 * (-pkin(5) + qJ(1)) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t55 - mrSges(1,2) + t65 * t57 - t64 * t47 + t63 * (t57 * pkin(2) + t55) + t62 * t37 + t61 * t36) * g(2) + (-m(3) * t58 - mrSges(1,1) + t64 * t57 + t65 * t47 + t63 * (t47 * pkin(2) + t58) - t61 * t37 + t62 * t36) * g(1);
U = t1;
