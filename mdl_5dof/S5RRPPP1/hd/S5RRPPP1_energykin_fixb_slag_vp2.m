% Calculate kinetic energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:17
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.32s
% Computational Cost: add. (362->88), mult. (976->121), div. (0->0), fcn. (678->6), ass. (0->32)
t99 = pkin(7) * mrSges(3,3);
t82 = sin(pkin(8));
t84 = cos(pkin(5));
t98 = t82 * t84;
t97 = pkin(3) + qJ(5);
t96 = cos(pkin(8));
t86 = cos(qJ(2));
t95 = qJD(1) * t86;
t83 = sin(pkin(5));
t94 = qJD(2) * t83;
t85 = sin(qJ(2));
t93 = t85 * qJD(1);
t76 = pkin(7) * t95 + (t84 * t95 + t94) * qJ(3);
t77 = qJD(2) * pkin(2) + (-qJ(3) * t84 - pkin(7)) * t93;
t78 = (-qJ(3) * t83 * t85 - pkin(2) * t86 - pkin(1)) * qJD(1);
t67 = t83 * t82 * t78 + t96 * t76 + t77 * t98;
t92 = t83 * t96;
t91 = t84 * t96;
t68 = -t83 * t77 + t84 * t78 + qJD(3);
t79 = t84 * qJD(2) - t83 * t95;
t65 = -t79 * qJ(4) - t67;
t70 = t82 * t94 + (t96 * t85 + t86 * t98) * qJD(1);
t90 = -t70 * qJ(4) + t68;
t66 = -t82 * t76 + t77 * t91 + t78 * t92;
t89 = qJD(4) - t66;
t69 = -qJD(2) * t92 + t82 * t93 - t91 * t95;
t64 = -t79 * pkin(3) + t89;
t63 = t69 * pkin(3) + t90;
t62 = -t69 * pkin(4) + qJD(5) - t65;
t61 = t97 * t69 + t90;
t60 = t70 * pkin(4) - t97 * t79 + t89;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + t64 * mrSges(5,2) + t62 * mrSges(6,2) - t65 * mrSges(5,3) - t60 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t79) * t79 + (t64 * mrSges(5,1) + t60 * mrSges(6,1) + t68 * mrSges(4,2) - t61 * mrSges(6,2) - t66 * mrSges(4,3) - t63 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t70 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t79) * t70 + (t68 * mrSges(4,1) + t65 * mrSges(5,1) - t62 * mrSges(6,1) - t63 * mrSges(5,2) - t67 * mrSges(4,3) + t61 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t69 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t79 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t70) * t69 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t85 ^ 2 + t86 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t99 + Ifges(3,2) / 0.2e1) * t86) * t86 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t86 + (t99 + Ifges(3,1) / 0.2e1) * t85) * t85) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t86 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t85) * qJD(2)) * qJD(1);
T = t1;
