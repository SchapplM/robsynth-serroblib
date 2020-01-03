% Calculate kinetic energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:33
% EndTime: 2019-12-31 17:47:33
% DurationCPUTime: 0.31s
% Computational Cost: add. (186->72), mult. (455->113), div. (0->0), fcn. (243->6), ass. (0->33)
t80 = sin(pkin(8));
t83 = cos(pkin(7));
t93 = t80 * t83;
t86 = qJD(1) ^ 2;
t92 = t86 * qJ(2) ^ 2;
t81 = sin(pkin(7));
t88 = -qJ(3) * t81 - pkin(1);
t64 = qJD(2) + ((-pkin(2) - qJ(4)) * t83 + t88) * qJD(1);
t89 = qJ(2) * qJD(1);
t72 = t81 * t89 + qJD(3);
t91 = qJD(1) * t81;
t69 = pkin(3) * t91 + t72;
t82 = cos(pkin(8));
t61 = t82 * t64 + t80 * t69;
t90 = qJD(1) * t83;
t70 = pkin(3) * t90 + t83 * t89 + qJD(4);
t60 = -t80 * t64 + t82 * t69;
t85 = cos(qJ(5));
t84 = sin(qJ(5));
t79 = t83 ^ 2;
t78 = t81 ^ 2;
t77 = -qJD(1) * pkin(1) + qJD(2);
t74 = t79 * t92;
t71 = t82 * t90 + qJD(5);
t68 = qJD(2) + (-pkin(2) * t83 + t88) * qJD(1);
t66 = (t81 * t84 - t85 * t93) * qJD(1);
t65 = (t81 * t85 + t84 * t93) * qJD(1);
t62 = (pkin(4) * t82 + pkin(6) * t80) * t90 + t70;
t59 = pkin(6) * t91 + t61;
t58 = -pkin(4) * t91 - t60;
t57 = t85 * t59 + t84 * t62;
t56 = -t84 * t59 + t85 * t62;
t1 = m(3) * (t77 ^ 2 + t78 * t92 + t74) / 0.2e1 + m(4) * (t68 ^ 2 + t72 ^ 2 + t74) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t71 + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,6) * t71 + Ifges(6,2) * t65 / 0.2e1) * t65 + ((t68 * mrSges(4,2) + t70 * (mrSges(5,1) * t82 - mrSges(5,2) * t80) - t77 * mrSges(3,1) + (t82 ^ 2 * Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + (Ifges(5,4) * t82 + Ifges(5,1) * t80 / 0.2e1) * t80) * t90 + (t60 * t80 - t61 * t82) * mrSges(5,3)) * t83 + (t72 * mrSges(4,1) + t60 * mrSges(5,1) - t68 * mrSges(4,3) - t61 * mrSges(5,2) + t77 * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t91 + (-Ifges(5,5) * t80 - Ifges(5,6) * t82 + Ifges(3,4) + Ifges(4,6)) * t90) * t81) * qJD(1) + (Ifges(2,3) / 0.2e1 + (mrSges(4,1) * t79 + (t78 + t79) * mrSges(3,3)) * qJ(2)) * t86;
T = t1;
