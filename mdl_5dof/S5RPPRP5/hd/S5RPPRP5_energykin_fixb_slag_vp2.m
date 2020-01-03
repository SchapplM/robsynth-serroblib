% Calculate kinetic energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.21s
% Computational Cost: add. (154->66), mult. (383->94), div. (0->0), fcn. (201->4), ass. (0->25)
t80 = qJD(1) ^ 2;
t84 = t80 * qJ(2) ^ 2;
t75 = sin(pkin(7));
t82 = t75 * qJD(1);
t66 = qJ(2) * t82 + qJD(3);
t64 = -pkin(6) * t82 + t66;
t76 = cos(pkin(7));
t83 = qJD(1) * t76;
t65 = (-pkin(6) + qJ(2)) * t83;
t78 = sin(qJ(4));
t79 = cos(qJ(4));
t57 = t78 * t64 + t79 * t65;
t72 = -qJD(1) * pkin(1) + qJD(2);
t60 = -pkin(2) * t83 - qJ(3) * t82 + t72;
t58 = pkin(3) * t83 - t60;
t56 = t79 * t64 - t78 * t65;
t74 = t76 ^ 2;
t73 = t75 ^ 2;
t69 = t74 * t84;
t63 = (t75 * t79 - t76 * t78) * qJD(1);
t62 = (t75 * t78 + t76 * t79) * qJD(1);
t55 = qJD(4) * qJ(5) + t57;
t54 = -qJD(4) * pkin(4) + qJD(5) - t56;
t53 = t62 * pkin(4) - t63 * qJ(5) + t58;
t1 = m(3) * (t72 ^ 2 + t73 * t84 + t69) / 0.2e1 + m(4) * (t60 ^ 2 + t66 ^ 2 + t69) / 0.2e1 + m(6) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t58 * mrSges(5,2) + t54 * mrSges(6,2) - t56 * mrSges(5,3) - t53 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t63) * t63 + (t58 * mrSges(5,1) + t53 * mrSges(6,1) - t55 * mrSges(6,2) - t57 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t62 + (-Ifges(5,4) + Ifges(6,5)) * t63) * t62 + (t56 * mrSges(5,1) - t54 * mrSges(6,1) - t57 * mrSges(5,2) + t55 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(4) + (Ifges(6,4) + Ifges(5,5)) * t63 + (-Ifges(5,6) + Ifges(6,6)) * t62) * qJD(4) + ((-t72 * mrSges(3,1) - t60 * mrSges(4,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t83) * t76 + (t72 * mrSges(3,2) + t66 * mrSges(4,2) - t60 * mrSges(4,3) + ((Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t75 + (Ifges(3,4) - Ifges(4,5)) * t76) * qJD(1)) * t75) * qJD(1) + (Ifges(2,3) / 0.2e1 + (mrSges(4,2) * t74 + (t73 + t74) * mrSges(3,3)) * qJ(2)) * t80;
T = t1;
