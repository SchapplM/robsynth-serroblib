% Calculate kinetic energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:55
% EndTime: 2019-12-31 18:03:55
% DurationCPUTime: 0.27s
% Computational Cost: add. (214->70), mult. (533->111), div. (0->0), fcn. (329->6), ass. (0->32)
t95 = qJD(1) ^ 2;
t99 = t95 * qJ(2) ^ 2;
t88 = sin(pkin(8));
t97 = t88 * qJD(1);
t78 = qJ(2) * t97 + qJD(3);
t76 = -pkin(6) * t97 + t78;
t89 = cos(pkin(8));
t98 = qJD(1) * t89;
t77 = (-pkin(6) + qJ(2)) * t98;
t92 = sin(qJ(4));
t94 = cos(qJ(4));
t68 = t92 * t76 + t94 * t77;
t84 = -qJD(1) * pkin(1) + qJD(2);
t67 = t94 * t76 - t92 * t77;
t72 = -pkin(2) * t98 - qJ(3) * t97 + t84;
t69 = pkin(3) * t98 - t72;
t93 = cos(qJ(5));
t91 = sin(qJ(5));
t87 = qJD(4) + qJD(5);
t86 = t89 ^ 2;
t85 = t88 ^ 2;
t81 = t86 * t99;
t75 = (t88 * t94 - t89 * t92) * qJD(1);
t74 = (-t88 * t92 - t89 * t94) * qJD(1);
t66 = t91 * t74 + t93 * t75;
t65 = t93 * t74 - t91 * t75;
t64 = -t74 * pkin(4) + t69;
t63 = t74 * pkin(7) + t68;
t62 = qJD(4) * pkin(4) - t75 * pkin(7) + t67;
t61 = t91 * t62 + t93 * t63;
t60 = t93 * t62 - t91 * t63;
t1 = m(3) * (t84 ^ 2 + t85 * t99 + t81) / 0.2e1 + m(4) * (t72 ^ 2 + t78 ^ 2 + t81) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t64 ^ 2) / 0.2e1 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t87 / 0.2e1) * t87 + (t69 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t75 / 0.2e1) * t75 + (-t69 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t75 + Ifges(5,2) * t74 / 0.2e1) * t74 + (t64 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t87 + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t64 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,6) * t87 + Ifges(6,2) * t65 / 0.2e1) * t65 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,5) * t75 + Ifges(5,6) * t74 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + ((-t84 * mrSges(3,1) - t72 * mrSges(4,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t98) * t89 + (t84 * mrSges(3,2) + t78 * mrSges(4,2) - t72 * mrSges(4,3) + ((Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t88 + (Ifges(3,4) - Ifges(4,5)) * t89) * qJD(1)) * t88) * qJD(1) + (Ifges(2,3) / 0.2e1 + (mrSges(4,2) * t86 + (t85 + t86) * mrSges(3,3)) * qJ(2)) * t95;
T = t1;
