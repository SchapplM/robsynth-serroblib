% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RPRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->9), mult. (41->23), div. (0->0), fcn. (41->6), ass. (0->14)
t116 = sin(qJ(1));
t126 = qJD(1) * t116;
t119 = cos(qJ(1));
t125 = qJD(1) * t119;
t115 = sin(qJ(3));
t124 = qJD(3) * t115;
t118 = cos(qJ(3));
t123 = qJD(3) * t118;
t122 = qJD(3) * t119;
t121 = qJD(1) * t118 - qJD(4);
t117 = cos(qJ(4));
t120 = (qJD(4) * t118 - qJD(1)) * t117;
t114 = sin(qJ(4));
t1 = [0, 0, t125, -t115 * t126 + t118 * t122, t119 * t120 + (-t115 * t122 - t121 * t116) * t114; 0, 0, t126, t115 * t125 + t116 * t123, t116 * t120 + (-t116 * t124 + t121 * t119) * t114; 0, 0, 0, t124, t115 * qJD(4) * t117 + t114 * t123;];
JgD_rot  = t1;
