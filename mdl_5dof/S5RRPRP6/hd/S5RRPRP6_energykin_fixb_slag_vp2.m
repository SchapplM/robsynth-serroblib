% Calculate kinetic energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:42
% EndTime: 2019-12-31 19:56:42
% DurationCPUTime: 0.37s
% Computational Cost: add. (318->82), mult. (752->119), div. (0->0), fcn. (486->6), ass. (0->27)
t91 = qJD(1) * (pkin(6) + qJ(3));
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t72 = (t79 * t82 - t80 * t84) * qJD(1);
t89 = pkin(6) * mrSges(3,3);
t73 = (t79 * t84 + t80 * t82) * qJD(1);
t78 = qJD(3) + (-pkin(2) * t84 - pkin(1)) * qJD(1);
t62 = pkin(3) * t72 - pkin(7) * t73 + t78;
t76 = qJD(2) * pkin(2) - t82 * t91;
t77 = t84 * t91;
t67 = t79 * t76 + t80 * t77;
t65 = qJD(2) * pkin(7) + t67;
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t58 = t81 * t62 + t83 * t65;
t57 = t83 * t62 - t65 * t81;
t66 = t76 * t80 - t79 * t77;
t64 = -qJD(2) * pkin(3) - t66;
t71 = qJD(4) + t72;
t69 = qJD(2) * t81 + t73 * t83;
t68 = qJD(2) * t83 - t73 * t81;
t59 = -pkin(4) * t68 + qJD(5) + t64;
t56 = qJ(5) * t68 + t58;
t55 = pkin(4) * t71 - qJ(5) * t69 + t57;
t1 = m(4) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (t78 * mrSges(4,2) - t66 * mrSges(4,3) + Ifges(4,1) * t73 / 0.2e1) * t73 - (-t78 * mrSges(4,1) + t67 * mrSges(4,3) + Ifges(4,4) * t73 - Ifges(4,2) * t72 / 0.2e1) * t72 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + Ifges(4,5) * t73 - Ifges(4,6) * t72 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t82 + Ifges(3,6) * t84 + (-mrSges(3,1) * t82 - mrSges(3,2) * t84) * pkin(6)) * qJD(1)) * qJD(2) + (t57 * mrSges(5,1) + t55 * mrSges(6,1) - t58 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t71) * t71 + (t64 * mrSges(5,2) + t59 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t69 + (Ifges(5,5) + Ifges(6,5)) * t71) * t69 + (-t64 * mrSges(5,1) - t59 * mrSges(6,1) + t58 * mrSges(5,3) + t56 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t68 + (Ifges(5,6) + Ifges(6,6)) * t71 + (Ifges(5,4) + Ifges(6,4)) * t69) * t68 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t82 ^ 2 + t84 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t89 + Ifges(3,2) / 0.2e1) * t84) * t84 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t84 + (t89 + Ifges(3,1) / 0.2e1) * t82) * t82) * qJD(1) ^ 2;
T = t1;
