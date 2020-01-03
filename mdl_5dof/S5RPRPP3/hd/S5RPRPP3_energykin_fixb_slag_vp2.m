% Calculate kinetic energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:17
% DurationCPUTime: 0.30s
% Computational Cost: add. (195->73), mult. (496->95), div. (0->0), fcn. (290->4), ass. (0->27)
t71 = cos(pkin(7));
t84 = t71 ^ 2;
t83 = m(3) / 0.2e1;
t82 = cos(qJ(3));
t81 = pkin(3) + qJ(5);
t80 = pkin(6) + qJ(2);
t70 = sin(pkin(7));
t77 = t70 * qJD(1);
t63 = t80 * t77;
t78 = qJD(1) * t71;
t64 = t80 * t78;
t72 = sin(qJ(3));
t57 = -t72 * t63 + t82 * t64;
t56 = -t82 * t63 - t72 * t64;
t55 = -qJD(3) * qJ(4) - t57;
t76 = qJD(4) - t56;
t65 = qJD(2) + (-pkin(2) * t71 - pkin(1)) * qJD(1);
t62 = (t82 * t70 + t71 * t72) * qJD(1);
t75 = -qJ(4) * t62 + t65;
t67 = -qJD(1) * pkin(1) + qJD(2);
t61 = t72 * t77 - t82 * t78;
t54 = -qJD(3) * pkin(3) + t76;
t53 = pkin(3) * t61 + t75;
t52 = -pkin(4) * t61 + qJD(5) - t55;
t51 = t62 * pkin(4) - t81 * qJD(3) + t76;
t50 = t81 * t61 + t75;
t1 = t67 ^ 2 * t83 + m(4) * (t56 ^ 2 + t57 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + (t54 * mrSges(5,1) + t51 * mrSges(6,1) + t65 * mrSges(4,2) - t50 * mrSges(6,2) - t56 * mrSges(4,3) - t53 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t62) * t62 + (t65 * mrSges(4,1) + t55 * mrSges(5,1) - t52 * mrSges(6,1) - t53 * mrSges(5,2) - t57 * mrSges(4,3) + t50 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t61 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t62) * t61 + (t56 * mrSges(4,1) - t57 * mrSges(4,2) + t54 * mrSges(5,2) + t52 * mrSges(6,2) - t55 * mrSges(5,3) - t51 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t62 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t61) * qJD(3) + (t67 * (-mrSges(3,1) * t71 + mrSges(3,2) * t70) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t83 + mrSges(3,3)) * (t70 ^ 2 + t84) * qJ(2) + Ifges(3,2) * t84 / 0.2e1 + (Ifges(3,4) * t71 + Ifges(3,1) * t70 / 0.2e1) * t70) * qJD(1)) * qJD(1);
T = t1;
