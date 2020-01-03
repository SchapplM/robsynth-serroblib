% Calculate kinetic energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:27
% EndTime: 2019-12-31 16:52:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (175->56), mult. (467->94), div. (0->0), fcn. (304->6), ass. (0->26)
t72 = cos(pkin(7));
t83 = t72 ^ 2;
t82 = qJD(1) * (pkin(5) + qJ(2));
t81 = m(3) / 0.2e1;
t71 = sin(pkin(7));
t64 = t71 * t82;
t65 = t72 * t82;
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t57 = -t74 * t64 + t76 * t65;
t56 = -t76 * t64 - t65 * t74;
t66 = qJD(2) + (-pkin(2) * t72 - pkin(1)) * qJD(1);
t75 = cos(qJ(4));
t73 = sin(qJ(4));
t70 = qJD(3) + qJD(4);
t67 = -qJD(1) * pkin(1) + qJD(2);
t63 = (t71 * t76 + t72 * t74) * qJD(1);
t62 = (-t71 * t74 + t72 * t76) * qJD(1);
t58 = -pkin(3) * t62 + t66;
t55 = t62 * t73 + t63 * t75;
t54 = t62 * t75 - t63 * t73;
t53 = pkin(6) * t62 + t57;
t52 = qJD(3) * pkin(3) - pkin(6) * t63 + t56;
t51 = t52 * t73 + t53 * t75;
t50 = t52 * t75 - t53 * t73;
t1 = t67 ^ 2 * t81 + m(5) * (t50 ^ 2 + t51 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t66 ^ 2) / 0.2e1 + (t50 * mrSges(5,1) - t51 * mrSges(5,2) + Ifges(5,3) * t70 / 0.2e1) * t70 + (t66 * mrSges(4,2) - t56 * mrSges(4,3) + Ifges(4,1) * t63 / 0.2e1) * t63 + (-t66 * mrSges(4,1) + t57 * mrSges(4,3) + Ifges(4,4) * t63 + Ifges(4,2) * t62 / 0.2e1) * t62 + (t58 * mrSges(5,2) - t50 * mrSges(5,3) + Ifges(5,5) * t70 + Ifges(5,1) * t55 / 0.2e1) * t55 + (-t58 * mrSges(5,1) + t51 * mrSges(5,3) + Ifges(5,4) * t55 + Ifges(5,6) * t70 + Ifges(5,2) * t54 / 0.2e1) * t54 + (t56 * mrSges(4,1) - t57 * mrSges(4,2) + Ifges(4,5) * t63 + Ifges(4,6) * t62 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t67 * (-mrSges(3,1) * t72 + mrSges(3,2) * t71) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t81 + mrSges(3,3)) * (t71 ^ 2 + t83) * qJ(2) + Ifges(3,2) * t83 / 0.2e1 + (Ifges(3,4) * t72 + Ifges(3,1) * t71 / 0.2e1) * t71) * qJD(1)) * qJD(1);
T = t1;
