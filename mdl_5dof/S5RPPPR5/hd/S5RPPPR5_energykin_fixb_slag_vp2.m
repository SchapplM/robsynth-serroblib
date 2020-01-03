% Calculate kinetic energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:19
% EndTime: 2019-12-31 17:46:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (160->54), mult. (306->87), div. (0->0), fcn. (140->6), ass. (0->26)
t88 = m(3) / 0.2e1;
t72 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t79 = sin(pkin(7));
t81 = cos(pkin(7));
t86 = qJ(2) * qJD(1);
t68 = t79 * t72 + t81 * t86;
t66 = -qJD(1) * qJ(4) + t68;
t78 = sin(pkin(8));
t80 = cos(pkin(8));
t62 = t78 * qJD(3) + t80 * t66;
t87 = t80 * qJD(1);
t67 = t81 * t72 - t79 * t86;
t65 = qJD(1) * pkin(3) + qJD(4) - t67;
t84 = cos(qJ(5));
t83 = sin(qJ(5));
t77 = t80 * qJD(3);
t75 = -qJD(1) * pkin(1) + qJD(2);
t70 = (-t78 * t84 - t80 * t83) * qJD(1);
t69 = (t78 * t83 - t80 * t84) * qJD(1);
t63 = pkin(4) * t87 + t65;
t61 = -t78 * t66 + t77;
t60 = -pkin(6) * t87 + t62;
t59 = t77 + (pkin(6) * qJD(1) - t66) * t78;
t58 = t83 * t59 + t84 * t60;
t57 = t84 * t59 - t83 * t60;
t1 = m(4) * (qJD(3) ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t63 ^ 2) / 0.2e1 + t75 ^ 2 * t88 + (t63 * mrSges(6,2) - t57 * mrSges(6,3) + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t63 * mrSges(6,1) + t58 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,5) * t70 + Ifges(6,6) * t69 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t68 * mrSges(4,2) - t67 * mrSges(4,1) + t65 * (mrSges(5,1) * t80 - mrSges(5,2) * t78) - t75 * mrSges(3,1) + (t61 * t78 - t62 * t80) * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + (qJ(2) * t88 + mrSges(3,3)) * qJ(2) + t80 ^ 2 * Ifges(5,2) / 0.2e1 + (Ifges(5,4) * t80 + Ifges(5,1) * t78 / 0.2e1) * t78) * qJD(1)) * qJD(1);
T = t1;
