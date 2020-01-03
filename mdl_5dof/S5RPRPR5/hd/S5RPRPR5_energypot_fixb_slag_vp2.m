% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:19
% EndTime: 2020-01-03 11:41:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (152->49), mult. (187->38), div. (0->0), fcn. (168->10), ass. (0->22)
t52 = qJ(3) + pkin(9);
t48 = qJ(5) + t52;
t43 = sin(t48);
t44 = cos(t48);
t58 = cos(qJ(3));
t45 = t58 * pkin(3) + pkin(2);
t46 = sin(t52);
t47 = cos(t52);
t56 = sin(qJ(3));
t77 = mrSges(3,1) + m(4) * pkin(2) + t58 * mrSges(4,1) - t56 * mrSges(4,2) + m(5) * t45 + t47 * mrSges(5,1) - t46 * mrSges(5,2) + m(6) * (pkin(4) * t47 + t45) + t44 * mrSges(6,1) - t43 * mrSges(6,2);
t55 = -qJ(4) - pkin(6);
t76 = -mrSges(3,2) + m(4) * pkin(6) - m(5) * t55 - m(6) * (-pkin(7) + t55) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t75 = m(3) + m(4);
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t70 = m(5) + m(6) + t75;
t74 = t70 * pkin(1) + t76 * t53 + t77 * t54 + mrSges(2,1);
t72 = -pkin(3) * t56 - qJ(2);
t69 = -m(5) * t72 + t46 * mrSges(5,1) + t47 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) - m(6) * (-pkin(4) * t46 + t72) + t43 * mrSges(6,1) + t44 * mrSges(6,2) + t56 * mrSges(4,1) + t58 * mrSges(4,2) + t75 * qJ(2);
t59 = cos(qJ(1));
t57 = sin(qJ(1));
t1 = (t69 * t57 + t74 * t59 - mrSges(1,3)) * g(3) + (-t74 * t57 + t69 * t59 - mrSges(1,2)) * g(2) + (-mrSges(1,1) - mrSges(2,3) + (-m(2) - t70) * pkin(5) + t76 * t54 - t77 * t53) * g(1);
U = t1;
