% Calculate potential energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:12
% EndTime: 2019-12-31 18:56:12
% DurationCPUTime: 0.34s
% Computational Cost: add. (96->50), mult. (150->46), div. (0->0), fcn. (127->6), ass. (0->21)
t56 = cos(qJ(4));
t79 = m(5) * pkin(3) + m(6) * (pkin(4) * t56 + pkin(3)) + mrSges(4,1);
t78 = -m(5) * pkin(7) + m(6) * (-qJ(5) - pkin(7)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t77 = m(6) * pkin(4);
t76 = -mrSges(5,1) - mrSges(6,1);
t75 = mrSges(5,2) + mrSges(6,2);
t74 = -m(4) - m(5) - m(6);
t73 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t72 = t79 * t54 + t78 * t57 - mrSges(2,2) + mrSges(3,3);
t53 = sin(qJ(4));
t58 = cos(qJ(1));
t70 = t53 * t58;
t55 = sin(qJ(1));
t69 = t55 * t53;
t68 = t55 * t56;
t65 = t58 * t54;
t63 = t58 * pkin(1) + t55 * qJ(2);
t49 = t55 * pkin(1);
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) + t74 * (pkin(2) + pkin(5)) + (t75 * t53 + t76 * t56 - t79) * t57 + t78 * t54) * g(3) + (-t69 * t77 - m(3) * t49 - mrSges(1,2) + t74 * (t55 * pkin(6) + t49) + t76 * (-t56 * t65 + t69) - t75 * (t53 * t65 + t68) + t73 * t55 + ((m(3) - t74) * qJ(2) + t72) * t58) * g(2) + (-t70 * t77 - m(3) * t63 - mrSges(1,1) + t74 * (t58 * pkin(6) + t63) + t76 * (t54 * t68 + t70) - t75 * (-t54 * t69 + t56 * t58) + t73 * t58 - t72 * t55) * g(1);
U = t1;
