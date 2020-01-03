% Calculate kinetic energy for
% S5RRRPP6
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:55
% EndTime: 2019-12-31 20:59:55
% DurationCPUTime: 0.33s
% Computational Cost: add. (354->83), mult. (766->119), div. (0->0), fcn. (488->6), ass. (0->28)
t90 = pkin(6) * mrSges(3,3);
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t72 = (-pkin(2) * t84 - pkin(7) * t82 - pkin(1)) * qJD(1);
t87 = t84 * qJD(1);
t77 = pkin(6) * t87 + qJD(2) * pkin(7);
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t66 = t83 * t72 - t77 * t81;
t88 = t82 * qJD(1);
t74 = qJD(2) * t81 + t83 * t88;
t78 = qJD(3) - t87;
t61 = pkin(3) * t78 - qJ(4) * t74 + t66;
t67 = t81 * t72 + t83 * t77;
t73 = qJD(2) * t83 - t81 * t88;
t63 = qJ(4) * t73 + t67;
t80 = sin(pkin(8));
t89 = cos(pkin(8));
t58 = t80 * t61 + t89 * t63;
t76 = -qJD(2) * pkin(2) + pkin(6) * t88;
t57 = t89 * t61 - t80 * t63;
t68 = -pkin(3) * t73 + qJD(4) + t76;
t65 = t80 * t73 + t89 * t74;
t64 = -t89 * t73 + t74 * t80;
t59 = pkin(4) * t64 - qJ(5) * t65 + t68;
t56 = qJ(5) * t78 + t58;
t55 = -t78 * pkin(4) + qJD(5) - t57;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (t76 * mrSges(4,2) - t66 * mrSges(4,3) + Ifges(4,1) * t74 / 0.2e1) * t74 + (-t76 * mrSges(4,1) + t67 * mrSges(4,3) + Ifges(4,4) * t74 + Ifges(4,2) * t73 / 0.2e1) * t73 + (t68 * mrSges(5,2) + t55 * mrSges(6,2) - t57 * mrSges(5,3) - t59 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t65) * t65 + (t68 * mrSges(5,1) + t59 * mrSges(6,1) - t56 * mrSges(6,2) - t58 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t64 + (-Ifges(5,4) + Ifges(6,5)) * t65) * t64 + (t66 * mrSges(4,1) + t57 * mrSges(5,1) - t55 * mrSges(6,1) - t67 * mrSges(4,2) - t58 * mrSges(5,2) + t56 * mrSges(6,3) + Ifges(4,5) * t74 + Ifges(4,6) * t73 + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t78 + (Ifges(6,4) + Ifges(5,5)) * t65 + (-Ifges(5,6) + Ifges(6,6)) * t64) * t78 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t82 ^ 2 + t84 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t90 + Ifges(3,2) / 0.2e1) * t84) * t84 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t84 + (t90 + Ifges(3,1) / 0.2e1) * t82) * t82) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t84 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t82) * qJD(2)) * qJD(1);
T = t1;
