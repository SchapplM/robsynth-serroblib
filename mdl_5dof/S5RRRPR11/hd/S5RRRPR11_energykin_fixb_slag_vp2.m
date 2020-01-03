% Calculate kinetic energy for
% S5RRRPR11
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:21
% EndTime: 2019-12-31 21:32:21
% DurationCPUTime: 0.36s
% Computational Cost: add. (304->84), mult. (622->121), div. (0->0), fcn. (366->6), ass. (0->32)
t99 = -pkin(3) - pkin(4);
t98 = pkin(6) * mrSges(3,3);
t97 = cos(qJ(3));
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t74 = (-pkin(2) * t90 - pkin(7) * t88 - pkin(1)) * qJD(1);
t95 = t90 * qJD(1);
t80 = pkin(6) * t95 + qJD(2) * pkin(7);
t87 = sin(qJ(3));
t72 = t87 * t74 + t97 * t80;
t96 = t88 * qJD(1);
t83 = qJD(3) - t95;
t67 = t83 * qJ(4) + t72;
t79 = -qJD(2) * pkin(2) + pkin(6) * t96;
t71 = t74 * t97 - t87 * t80;
t94 = qJD(4) - t71;
t76 = t87 * qJD(2) + t96 * t97;
t93 = qJ(4) * t76 - t79;
t89 = cos(qJ(5));
t86 = sin(qJ(5));
t82 = qJD(5) - t83;
t75 = -qJD(2) * t97 + t87 * t96;
t70 = t75 * t86 + t76 * t89;
t69 = t75 * t89 - t76 * t86;
t68 = pkin(3) * t75 - t93;
t66 = -t83 * pkin(3) + t94;
t65 = t75 * t99 + t93;
t64 = pkin(8) * t75 + t67;
t63 = -t76 * pkin(8) + t83 * t99 + t94;
t62 = t63 * t86 + t64 * t89;
t61 = t63 * t89 - t64 * t86;
t1 = m(4) * (t71 ^ 2 + t72 ^ 2 + t79 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t65 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t65 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t82 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t71 * mrSges(4,1) - t66 * mrSges(5,1) - t72 * mrSges(4,2) + t67 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t83) * t83 + (t79 * mrSges(4,2) + t66 * mrSges(5,2) - t71 * mrSges(4,3) - t68 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t76 + (Ifges(5,4) + Ifges(4,5)) * t83) * t76 + (t79 * mrSges(4,1) + t68 * mrSges(5,1) - t67 * mrSges(5,2) - t72 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t75 + (-Ifges(4,6) + Ifges(5,6)) * t83 + (-Ifges(4,4) + Ifges(5,5)) * t76) * t75 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t88 ^ 2 + t90 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t98 + Ifges(3,2) / 0.2e1) * t90) * t90 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t90 + (t98 + Ifges(3,1) / 0.2e1) * t88) * t88) * qJD(1) + ((-mrSges(3,2) * pkin(6) + Ifges(3,6)) * t90 + (-mrSges(3,1) * pkin(6) + Ifges(3,5)) * t88) * qJD(2)) * qJD(1);
T = t1;
