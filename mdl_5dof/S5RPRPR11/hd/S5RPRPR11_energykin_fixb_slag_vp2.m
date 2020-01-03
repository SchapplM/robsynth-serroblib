% Calculate kinetic energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:06
% EndTime: 2019-12-31 18:27:06
% DurationCPUTime: 0.35s
% Computational Cost: add. (255->77), mult. (642->112), div. (0->0), fcn. (410->6), ass. (0->34)
t101 = cos(qJ(3));
t87 = sin(pkin(8));
t88 = cos(pkin(8));
t90 = sin(qJ(3));
t77 = (t101 * t87 + t88 * t90) * qJD(1);
t80 = -(pkin(2) * t88 + pkin(1)) * qJD(1) + qJD(2);
t105 = -t77 * qJ(4) + t80;
t104 = t88 ^ 2;
t103 = m(3) / 0.2e1;
t102 = -pkin(3) - pkin(4);
t100 = pkin(6) + qJ(2);
t96 = t87 * qJD(1);
t78 = t100 * t96;
t97 = qJD(1) * t88;
t79 = t100 * t97;
t72 = t101 * t79 - t90 * t78;
t70 = qJD(3) * qJ(4) + t72;
t71 = -t101 * t78 - t90 * t79;
t94 = qJD(4) - t71;
t91 = cos(qJ(5));
t89 = sin(qJ(5));
t85 = -qJD(3) + qJD(5);
t82 = -qJD(1) * pkin(1) + qJD(2);
t76 = -t101 * t97 + t90 * t96;
t69 = -qJD(3) * pkin(3) + t94;
t68 = t89 * t76 + t91 * t77;
t67 = t91 * t76 - t89 * t77;
t66 = t76 * pkin(3) + t105;
t65 = t76 * pkin(7) + t70;
t64 = -t77 * pkin(7) + t102 * qJD(3) + t94;
t63 = t102 * t76 - t105;
t62 = t89 * t64 + t91 * t65;
t61 = t91 * t64 - t89 * t65;
t1 = t82 ^ 2 * t103 + m(4) * (t71 ^ 2 + t72 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t85 / 0.2e1) * t85 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t85 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t85 + Ifges(6,2) * t67 / 0.2e1) * t67 + (t80 * mrSges(4,2) + t69 * mrSges(5,2) - t71 * mrSges(4,3) - t66 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t77) * t77 + (t80 * mrSges(4,1) + t66 * mrSges(5,1) - t70 * mrSges(5,2) - t72 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t76 + (-Ifges(4,4) + Ifges(5,5)) * t77) * t76 + (t71 * mrSges(4,1) - t69 * mrSges(5,1) - t72 * mrSges(4,2) + t70 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(5,4) + Ifges(4,5)) * t77 + (-Ifges(4,6) + Ifges(5,6)) * t76) * qJD(3) + (t82 * (-mrSges(3,1) * t88 + mrSges(3,2) * t87) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t103 + mrSges(3,3)) * (t87 ^ 2 + t104) * qJ(2) + Ifges(3,2) * t104 / 0.2e1 + (Ifges(3,4) * t88 + Ifges(3,1) * t87 / 0.2e1) * t87) * qJD(1)) * qJD(1);
T = t1;
