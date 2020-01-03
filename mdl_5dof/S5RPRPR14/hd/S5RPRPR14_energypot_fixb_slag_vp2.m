% Calculate potential energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:27
% EndTime: 2019-12-31 18:34:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (108->52), mult. (130->40), div. (0->0), fcn. (103->8), ass. (0->20)
t67 = m(5) + m(6);
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t66 = -m(6) * pkin(4) - t48 * mrSges(6,1) + t45 * mrSges(6,2) - mrSges(5,1);
t65 = m(3) + m(4);
t63 = -m(6) * pkin(7) - mrSges(6,3);
t43 = qJ(3) + pkin(8);
t37 = sin(t43);
t38 = cos(t43);
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t62 = -t46 * mrSges(4,1) - t49 * mrSges(4,2) - t38 * mrSges(5,2) + t66 * t37 + mrSges(2,2) - mrSges(3,3);
t61 = -m(4) * pkin(6) - t45 * mrSges(6,1) - t48 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + t67 * (-qJ(4) - pkin(6));
t60 = pkin(2) + pkin(5);
t59 = pkin(3) * t46;
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t58 = t50 * pkin(1) + t47 * qJ(2);
t55 = -qJ(2) - t59;
t1 = (-m(4) * t60 - t49 * mrSges(4,1) + t46 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t67 * (t49 * pkin(3) + t60) + t66 * t38 + (mrSges(5,2) + t63) * t37 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + (-m(5) * t55 - m(6) * (pkin(7) * t38 + t55) - t38 * mrSges(6,3) + t65 * qJ(2) - t62) * t50 + ((-t67 - t65) * pkin(1) + t61) * t47) * g(2) + (-mrSges(1,1) - t65 * t58 - t67 * (t47 * t59 + t58) + t61 * t50 + (-t63 * t38 + t62) * t47) * g(1);
U = t1;
