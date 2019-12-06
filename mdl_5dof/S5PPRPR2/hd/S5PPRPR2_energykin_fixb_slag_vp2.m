% Calculate kinetic energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:54
% EndTime: 2019-12-05 15:02:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (71->39), mult. (155->62), div. (0->0), fcn. (78->6), ass. (0->18)
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t63 = (t65 * t70 + t66 * t68) * qJD(1);
t60 = qJD(3) * qJ(4) + t63;
t76 = t60 ^ 2;
t62 = (-t65 * t68 + t66 * t70) * qJD(1);
t74 = qJD(4) - t62;
t73 = qJD(1) ^ 2;
t72 = qJD(2) ^ 2;
t69 = cos(qJ(5));
t67 = sin(qJ(5));
t59 = -qJD(3) * pkin(3) + t74;
t58 = (-pkin(3) - pkin(6)) * qJD(3) + t74;
t57 = t69 * qJD(2) + t67 * t58;
t56 = -t67 * qJD(2) + t69 * t58;
t1 = m(3) * (t72 + (t65 ^ 2 + t66 ^ 2) * t73) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t72) / 0.2e1 + m(2) * t73 / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t76) / 0.2e1 + m(5) * (t59 ^ 2 + t72 + t76) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t62 * mrSges(4,1) - t63 * mrSges(4,2) + t59 * mrSges(5,2) + t60 * mrSges(5,3) + (t60 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t69 + (t60 * mrSges(6,1) - t57 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t67 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) * t69 ^ 2 / 0.2e1 + (-Ifges(6,4) * t69 + Ifges(6,2) * t67 / 0.2e1) * t67) * qJD(3)) * qJD(3);
T = t1;
