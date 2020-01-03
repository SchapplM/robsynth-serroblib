% Calculate kinetic energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR15_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:34
% EndTime: 2019-12-31 18:36:34
% DurationCPUTime: 0.24s
% Computational Cost: add. (239->71), mult. (483->114), div. (0->0), fcn. (266->6), ass. (0->29)
t76 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t89 = t76 * mrSges(4,3);
t83 = sin(qJ(3));
t85 = cos(qJ(3));
t71 = (pkin(3) * t83 - qJ(4) * t85 + qJ(2)) * qJD(1);
t72 = qJD(3) * qJ(4) + t76 * t83;
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t63 = t80 * t71 + t81 * t72;
t88 = qJD(1) * t83;
t87 = qJD(1) * t85;
t62 = t81 * t71 - t72 * t80;
t70 = -qJD(3) * pkin(3) - t76 * t85 + qJD(4);
t86 = qJD(1) ^ 2;
t84 = cos(qJ(5));
t82 = sin(qJ(5));
t79 = t86 * qJ(2) ^ 2;
t78 = -qJD(1) * pkin(1) + qJD(2);
t77 = qJD(5) + t88;
t74 = qJD(3) * t80 + t81 * t87;
t73 = qJD(3) * t81 - t80 * t87;
t66 = -pkin(4) * t73 + t70;
t65 = t73 * t82 + t74 * t84;
t64 = t73 * t84 - t74 * t82;
t61 = pkin(7) * t73 + t63;
t60 = pkin(4) * t88 - pkin(7) * t74 + t62;
t59 = t60 * t82 + t61 * t84;
t58 = t60 * t84 - t61 * t82;
t1 = m(6) * (t58 ^ 2 + t59 ^ 2 + t66 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t70 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t79) / 0.2e1 + m(4) * (t79 + (t83 ^ 2 + t85 ^ 2) * t76 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t86 + (t58 * mrSges(6,1) - t59 * mrSges(6,2) + Ifges(6,3) * t77 / 0.2e1) * t77 + (t70 * mrSges(5,2) - t62 * mrSges(5,3) + Ifges(5,1) * t74 / 0.2e1) * t74 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t85 * mrSges(4,1) - t83 * mrSges(4,2)) * t76) * qJD(3) + (t66 * mrSges(6,2) - t58 * mrSges(6,3) + Ifges(6,5) * t77 + Ifges(6,1) * t65 / 0.2e1) * t65 + (t78 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (-t89 + Ifges(4,1) * qJD(1) / 0.2e1) * t85) * t85 + (t74 * Ifges(5,5) - Ifges(4,6) * qJD(3) - t63 * mrSges(5,2) + t62 * mrSges(5,1) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t85) * qJD(1) + (-t89 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t83) * t83) * qJD(1) + (Ifges(5,6) * t88 - t70 * mrSges(5,1) + t63 * mrSges(5,3) + Ifges(5,4) * t74 + Ifges(5,2) * t73 / 0.2e1) * t73 + (-t66 * mrSges(6,1) + t59 * mrSges(6,3) + Ifges(6,4) * t65 + Ifges(6,6) * t77 + Ifges(6,2) * t64 / 0.2e1) * t64;
T = t1;
