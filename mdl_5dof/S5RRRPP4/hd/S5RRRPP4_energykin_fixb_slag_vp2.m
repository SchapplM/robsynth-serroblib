% Calculate kinetic energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:48
% EndTime: 2019-12-31 20:54:49
% DurationCPUTime: 0.36s
% Computational Cost: add. (354->82), mult. (844->118), div. (0->0), fcn. (566->6), ass. (0->27)
t90 = qJD(1) * (-pkin(7) - pkin(6));
t88 = pkin(6) * mrSges(3,3);
t82 = sin(qJ(2));
t76 = qJD(2) * pkin(2) + t82 * t90;
t84 = cos(qJ(2));
t77 = t84 * t90;
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t67 = t83 * t76 + t77 * t81;
t74 = (t81 * t84 + t82 * t83) * qJD(1);
t79 = qJD(2) + qJD(3);
t62 = pkin(3) * t79 - qJ(4) * t74 + t67;
t68 = t81 * t76 - t83 * t77;
t73 = (-t81 * t82 + t83 * t84) * qJD(1);
t64 = qJ(4) * t73 + t68;
t80 = sin(pkin(8));
t87 = cos(pkin(8));
t59 = t80 * t62 + t87 * t64;
t78 = (-pkin(2) * t84 - pkin(1)) * qJD(1);
t58 = t87 * t62 - t80 * t64;
t69 = -pkin(3) * t73 + qJD(4) + t78;
t66 = t80 * t73 + t87 * t74;
t65 = -t87 * t73 + t74 * t80;
t60 = pkin(4) * t65 - qJ(5) * t66 + t69;
t57 = qJ(5) * t79 + t59;
t56 = -t79 * pkin(4) + qJD(5) - t58;
t1 = m(5) * (t58 ^ 2 + t59 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t67 ^ 2 + t68 ^ 2 + t78 ^ 2) / 0.2e1 + (t78 * mrSges(4,2) - t67 * mrSges(4,3) + Ifges(4,1) * t74 / 0.2e1) * t74 + (-t78 * mrSges(4,1) + t68 * mrSges(4,3) + Ifges(4,4) * t74 + Ifges(4,2) * t73 / 0.2e1) * t73 + (t69 * mrSges(5,2) + t56 * mrSges(6,2) - t58 * mrSges(5,3) - t60 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t66) * t66 + (t69 * mrSges(5,1) + t60 * mrSges(6,1) - t57 * mrSges(6,2) - t59 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t65 + (-Ifges(5,4) + Ifges(6,5)) * t66) * t65 + (t67 * mrSges(4,1) + t58 * mrSges(5,1) - t56 * mrSges(6,1) - t68 * mrSges(4,2) - t59 * mrSges(5,2) + t57 * mrSges(6,3) + Ifges(4,5) * t74 + Ifges(4,6) * t73 + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t79 + (Ifges(6,4) + Ifges(5,5)) * t66 + (-Ifges(5,6) + Ifges(6,6)) * t65) * t79 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t82 ^ 2 + t84 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t88) * t84) * t84 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t84 + (Ifges(3,1) / 0.2e1 + t88) * t82) * t82) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t84 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t82) * qJD(2)) * qJD(1);
T = t1;
