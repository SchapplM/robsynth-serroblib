% Calculate kinetic energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:28
% EndTime: 2019-12-31 18:58:28
% DurationCPUTime: 0.20s
% Computational Cost: add. (183->67), mult. (347->98), div. (0->0), fcn. (160->4), ass. (0->23)
t76 = cos(qJ(4));
t64 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t75 = t64 * mrSges(4,3);
t74 = qJD(1) / 0.2e1;
t70 = sin(qJ(3));
t71 = cos(qJ(3));
t58 = (pkin(3) * t70 - pkin(7) * t71 + qJ(2)) * qJD(1);
t59 = qJD(3) * pkin(7) + t70 * t64;
t69 = sin(qJ(4));
t55 = t69 * t58 + t76 * t59;
t73 = qJD(1) * t71;
t60 = -qJD(3) * pkin(3) - t71 * t64;
t54 = t76 * t58 - t69 * t59;
t72 = qJD(1) ^ 2;
t68 = t72 * qJ(2) ^ 2;
t66 = -qJD(1) * pkin(1) + qJD(2);
t65 = t70 * qJD(1) + qJD(4);
t62 = t69 * qJD(3) + t76 * t73;
t61 = -t76 * qJD(3) + t69 * t73;
t53 = t61 * pkin(4) - t62 * qJ(5) + t60;
t52 = t65 * qJ(5) + t55;
t51 = -t65 * pkin(4) + qJD(5) - t54;
t1 = m(3) * (t66 ^ 2 + t68) / 0.2e1 + m(4) * (t68 + (t70 ^ 2 + t71 ^ 2) * t64 ^ 2) / 0.2e1 + m(5) * (t54 ^ 2 + t55 ^ 2 + t60 ^ 2) / 0.2e1 + m(6) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t72 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t71 * mrSges(4,1) - t70 * mrSges(4,2)) * t64) * qJD(3) + (t66 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t74 - t75) * t71) * t71 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t71) * qJD(1) + (Ifges(4,2) * t74 - t75) * t70) * t70) * qJD(1) + (t54 * mrSges(5,1) - t51 * mrSges(6,1) - t55 * mrSges(5,2) + t52 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t65) * t65 + (t60 * mrSges(5,2) + t51 * mrSges(6,2) - t54 * mrSges(5,3) - t53 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t62 + (Ifges(6,4) + Ifges(5,5)) * t65) * t62 + (t60 * mrSges(5,1) + t53 * mrSges(6,1) - t52 * mrSges(6,2) - t55 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t61 + (-Ifges(5,6) + Ifges(6,6)) * t65 + (-Ifges(5,4) + Ifges(6,5)) * t62) * t61;
T = t1;
