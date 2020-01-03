% Calculate kinetic energy for
% S5RRRPR8
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:02
% EndTime: 2019-12-31 21:19:02
% DurationCPUTime: 0.36s
% Computational Cost: add. (302->83), mult. (648->120), div. (0->0), fcn. (396->6), ass. (0->33)
t95 = pkin(3) + pkin(8);
t94 = -pkin(7) - pkin(6);
t93 = pkin(6) * mrSges(3,3);
t92 = cos(qJ(3));
t83 = sin(qJ(2));
t91 = t83 * qJD(1);
t76 = qJD(2) * pkin(2) + t91 * t94;
t85 = cos(qJ(2));
t90 = t85 * qJD(1);
t77 = t94 * t90;
t82 = sin(qJ(3));
t67 = t82 * t76 - t92 * t77;
t80 = qJD(2) + qJD(3);
t65 = -qJ(4) * t80 - t67;
t66 = t76 * t92 + t82 * t77;
t78 = (-pkin(2) * t85 - pkin(1)) * qJD(1);
t89 = qJD(4) - t66;
t73 = (t82 * t85 + t83 * t92) * qJD(1);
t88 = -qJ(4) * t73 + t78;
t84 = cos(qJ(5));
t81 = sin(qJ(5));
t72 = t82 * t91 - t90 * t92;
t71 = qJD(5) + t73;
t69 = t72 * t81 + t80 * t84;
t68 = t72 * t84 - t80 * t81;
t64 = -t80 * pkin(3) + t89;
t63 = pkin(3) * t72 + t88;
t62 = -pkin(4) * t72 - t65;
t61 = t72 * t95 + t88;
t60 = t73 * pkin(4) - t80 * t95 + t89;
t59 = t60 * t81 + t61 * t84;
t58 = t60 * t84 - t61 * t81;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t58 ^ 2 + t59 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + (t58 * mrSges(6,1) - t59 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (t62 * mrSges(6,2) - t58 * mrSges(6,3) + Ifges(6,5) * t71 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t62 * mrSges(6,1) + t59 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t71 + Ifges(6,2) * t68 / 0.2e1) * t68 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + t64 * mrSges(5,2) - t65 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t80) * t80 + (t64 * mrSges(5,1) + t78 * mrSges(4,2) - t66 * mrSges(4,3) - t63 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t73 + (-Ifges(5,4) + Ifges(4,5)) * t80) * t73 + (t78 * mrSges(4,1) + t65 * mrSges(5,1) - t63 * mrSges(5,2) - t67 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t72 + (Ifges(5,5) - Ifges(4,6)) * t80 + (-Ifges(4,4) - Ifges(5,6)) * t73) * t72 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t83 ^ 2 + t85 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t93 + Ifges(3,2) / 0.2e1) * t85) * t85 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t85 + (t93 + Ifges(3,1) / 0.2e1) * t83) * t83) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t85 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t83) * qJD(2)) * qJD(1);
T = t1;
