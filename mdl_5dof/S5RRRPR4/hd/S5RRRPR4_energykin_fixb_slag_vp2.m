% Calculate kinetic energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:10:56
% EndTime: 2019-12-31 21:10:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (234->70), mult. (329->108), div. (0->0), fcn. (140->6), ass. (0->28)
t93 = pkin(3) + pkin(4);
t77 = qJD(1) + qJD(2);
t81 = sin(qJ(2));
t90 = pkin(1) * qJD(1);
t74 = t77 * pkin(7) + t81 * t90;
t92 = t74 * mrSges(4,3);
t83 = cos(qJ(3));
t91 = t77 * t83;
t70 = qJD(3) * qJ(4) + t83 * t74;
t80 = sin(qJ(3));
t89 = t80 * t74 + qJD(4);
t84 = cos(qJ(2));
t88 = t84 * t90;
t87 = qJ(4) * t80 + pkin(2);
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t76 = -qJD(3) + qJD(5);
t75 = -t77 * pkin(2) - t88;
t69 = (-t79 * t83 + t80 * t82) * t77;
t68 = (-t79 * t80 - t82 * t83) * t77;
t67 = -qJD(3) * pkin(3) + t89;
t66 = -pkin(8) * t91 + t70;
t65 = -t88 + (-pkin(3) * t83 - t87) * t77;
t64 = -t80 * t77 * pkin(8) - t93 * qJD(3) + t89;
t63 = t88 + (t93 * t83 + t87) * t77;
t62 = t79 * t64 + t82 * t66;
t61 = t82 * t64 - t79 * t66;
t1 = m(4) * (t75 ^ 2 + (t80 ^ 2 + t83 ^ 2) * t74 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t81 ^ 2 + t84 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t76 / 0.2e1) * t76 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t76 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t76 + Ifges(6,2) * t68 / 0.2e1) * t68 + (-t67 * mrSges(5,1) + t70 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (-t80 * mrSges(4,1) - t83 * mrSges(4,2)) * t74) * qJD(3) + (Ifges(3,3) * t77 / 0.2e1 + (mrSges(3,1) * t84 - mrSges(3,2) * t81) * t90 + (-t75 * mrSges(4,1) - t65 * mrSges(5,1) + t70 * mrSges(5,2) + (t92 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t77) * t83 + (Ifges(4,6) - Ifges(5,6)) * qJD(3)) * t83 + (t67 * mrSges(5,2) - t65 * mrSges(5,3) + t75 * mrSges(4,2) + (t92 + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t77) * t80 + (Ifges(4,4) - Ifges(5,5)) * t91 + (Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t80) * t77;
T = t1;
