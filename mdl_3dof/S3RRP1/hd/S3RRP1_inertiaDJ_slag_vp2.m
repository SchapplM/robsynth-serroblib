% Calculate time derivative of joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% MqD [3x3]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaDJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_inertiaDJ_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaDJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_inertiaDJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_inertiaDJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_inertiaDJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:51
% EndTime: 2019-03-08 18:06:51
% DurationCPUTime: 0.05s
% Computational Cost: add. (20->13), mult. (62->23), div. (0->0), fcn. (16->2), ass. (0->10)
t9 = pkin(1) * qJD(2);
t4 = sin(qJ(2));
t8 = t4 * t9;
t5 = cos(qJ(2));
t7 = t5 * t9;
t1 = qJD(3) + t7;
t6 = -mrSges(3,2) * t7 + t1 * mrSges(4,3) + (-mrSges(3,1) - mrSges(4,1)) * t8;
t3 = qJD(3) * mrSges(4,3);
t2 = t4 * pkin(1) + qJ(3);
t10 = [0.2e1 * m(4) * (t2 * t1 + (-t5 * pkin(1) - pkin(2)) * t8) + 0.2e1 * t6; t3 + m(4) * (-pkin(2) * t8 + qJ(3) * t1 + qJD(3) * t2) + t6; 0.2e1 * m(4) * qJ(3) * qJD(3) + 0.2e1 * t3; m(4) * t8; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t10(1) t10(2) t10(4); t10(2) t10(3) t10(5); t10(4) t10(5) t10(6);];
Mq  = res;
