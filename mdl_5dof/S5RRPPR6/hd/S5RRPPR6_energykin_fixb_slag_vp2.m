% Calculate kinetic energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:31
% EndTime: 2019-12-31 19:31:31
% DurationCPUTime: 0.42s
% Computational Cost: add. (416->86), mult. (1002->134), div. (0->0), fcn. (690->8), ass. (0->35)
t104 = pkin(6) * mrSges(3,3);
t103 = pkin(6) + qJ(3);
t97 = cos(qJ(2));
t100 = qJD(1) * t97;
t95 = sin(qJ(2));
t101 = qJD(1) * t95;
t102 = cos(pkin(8));
t92 = sin(pkin(8));
t83 = -t102 * t100 + t101 * t92;
t84 = (t102 * t95 + t92 * t97) * qJD(1);
t89 = qJD(3) + (-pkin(2) * t97 - pkin(1)) * qJD(1);
t73 = pkin(3) * t83 - qJ(4) * t84 + t89;
t87 = qJD(2) * pkin(2) - t101 * t103;
t88 = t103 * t100;
t78 = t102 * t88 + t92 * t87;
t76 = qJD(2) * qJ(4) + t78;
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t67 = t91 * t73 + t93 * t76;
t66 = t93 * t73 - t76 * t91;
t77 = t102 * t87 - t92 * t88;
t75 = -qJD(2) * pkin(3) + qJD(4) - t77;
t96 = cos(qJ(5));
t94 = sin(qJ(5));
t82 = qJD(5) + t83;
t80 = qJD(2) * t91 + t84 * t93;
t79 = qJD(2) * t93 - t84 * t91;
t70 = t79 * t94 + t80 * t96;
t69 = t79 * t96 - t80 * t94;
t68 = -t79 * pkin(4) + t75;
t65 = pkin(7) * t79 + t67;
t64 = pkin(4) * t83 - pkin(7) * t80 + t66;
t63 = t64 * t94 + t65 * t96;
t62 = t64 * t96 - t65 * t94;
t1 = m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t89 ^ 2) / 0.2e1 + (t89 * mrSges(4,2) - t77 * mrSges(4,3) + Ifges(4,1) * t84 / 0.2e1) * t84 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t75 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (-t75 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,2) * t79 / 0.2e1) * t79 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t82 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t77 * mrSges(4,1) - t78 * mrSges(4,2) + Ifges(4,5) * t84 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t95 + Ifges(3,6) * t97 + (-mrSges(3,1) * t95 - mrSges(3,2) * t97) * pkin(6)) * qJD(1)) * qJD(2) + (t89 * mrSges(4,1) + t66 * mrSges(5,1) - t67 * mrSges(5,2) - t78 * mrSges(4,3) - Ifges(4,4) * t84 + Ifges(5,5) * t80 - Ifges(4,6) * qJD(2) + Ifges(5,6) * t79 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t83) * t83 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t95 ^ 2 + t97 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t104 + Ifges(3,2) / 0.2e1) * t97) * t97 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t97 + (t104 + Ifges(3,1) / 0.2e1) * t95) * t95) * qJD(1) ^ 2;
T = t1;
