% Calculate kinetic energy for
% S3PPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3PPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_energykin_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR2_energykin_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_energykin_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR2_energykin_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPR2_energykin_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPR2_energykin_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:00
% EndTime: 2018-11-14 10:11:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->9), mult. (26->19), div. (0->0), fcn. (8->2), ass. (0->6)
t25 = qJD(1) ^ 2;
t24 = cos(qJ(3));
t23 = sin(qJ(3));
t22 = t24 * qJD(1) + t23 * qJD(2);
t21 = -t23 * qJD(1) + t24 * qJD(2);
t1 = m(2) * t25 / 0.2e1 + m(3) * (qJD(2) ^ 2 + t25) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2) / 0.2e1 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3);
T  = t1;
