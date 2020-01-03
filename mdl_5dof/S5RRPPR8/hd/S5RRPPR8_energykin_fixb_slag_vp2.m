% Calculate kinetic energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:03
% DurationCPUTime: 0.36s
% Computational Cost: add. (264->84), mult. (594->124), div. (0->0), fcn. (336->6), ass. (0->30)
t100 = pkin(6) * mrSges(3,3);
t91 = sin(qJ(2));
t96 = qJ(4) * qJD(1);
t99 = qJD(1) * t91;
t97 = pkin(6) * t99 + qJD(3);
t72 = -t91 * t96 + (-pkin(2) - pkin(3)) * qJD(2) + t97;
t93 = cos(qJ(2));
t98 = qJD(1) * t93;
t79 = pkin(6) * t98 + qJD(2) * qJ(3);
t76 = -t93 * t96 + t79;
t87 = sin(pkin(8));
t88 = cos(pkin(8));
t65 = t87 * t72 + t88 * t76;
t77 = -qJD(1) * pkin(1) - pkin(2) * t98 - qJ(3) * t99;
t64 = t88 * t72 - t76 * t87;
t71 = pkin(3) * t98 + qJD(4) - t77;
t92 = cos(qJ(5));
t90 = sin(qJ(5));
t85 = -qJD(2) + qJD(5);
t78 = -qJD(2) * pkin(2) + t97;
t75 = (-t87 * t93 + t88 * t91) * qJD(1);
t74 = (-t87 * t91 - t88 * t93) * qJD(1);
t68 = t74 * t90 + t75 * t92;
t67 = t74 * t92 - t75 * t90;
t66 = -pkin(4) * t74 + t71;
t63 = pkin(7) * t74 + t65;
t62 = -qJD(2) * pkin(4) - pkin(7) * t75 + t64;
t61 = t62 * t90 + t63 * t92;
t60 = t62 * t92 - t63 * t90;
t1 = m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t66 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t71 ^ 2) / 0.2e1 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t85 / 0.2e1) * t85 + (t71 * mrSges(5,2) - t64 * mrSges(5,3) + Ifges(5,1) * t75 / 0.2e1) * t75 + (-t71 * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t75 + Ifges(5,2) * t74 / 0.2e1) * t74 + (t66 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t85 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t66 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t85 + Ifges(6,2) * t67 / 0.2e1) * t67 + (-t78 * mrSges(4,1) - t64 * mrSges(5,1) + t65 * mrSges(5,2) + t79 * mrSges(4,3) - Ifges(5,5) * t75 - Ifges(5,6) * t74 + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t77 * mrSges(4,1) + t79 * mrSges(4,2)) * t93 + (t78 * mrSges(4,2) - t77 * mrSges(4,3)) * t91 + ((-mrSges(3,2) * pkin(6) + Ifges(3,6) - Ifges(4,6)) * t93 + (-mrSges(3,1) * pkin(6) + Ifges(4,4) + Ifges(3,5)) * t91) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t91 ^ 2 + t93 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t100 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t93) * t93 + (-pkin(1) * mrSges(3,2) + (t100 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t91 + (Ifges(3,4) - Ifges(4,5)) * t93) * t91) * qJD(1)) * qJD(1);
T = t1;
