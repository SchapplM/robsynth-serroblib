% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_4_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_4_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_4_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (12->10), mult. (41->23), div. (0->0), fcn. (41->6), ass. (0->14)
t115 = sin(qJ(1));
t125 = qJD(1) * t115;
t118 = cos(qJ(1));
t124 = qJD(1) * t118;
t114 = sin(qJ(2));
t123 = qJD(2) * t114;
t117 = cos(qJ(2));
t122 = qJD(2) * t117;
t121 = qJD(2) * t118;
t120 = qJD(1) * t117 + qJD(3);
t116 = cos(qJ(3));
t119 = (-qJD(3) * t117 - qJD(1)) * t116;
t113 = sin(qJ(3));
t1 = [0, t124, t114 * t125 - t117 * t121, t118 * t119 + (t114 * t121 + t120 * t115) * t113, 0, 0, 0; 0, t125, -t114 * t124 - t115 * t122, t115 * t119 + (t115 * t123 - t120 * t118) * t113, 0, 0, 0; 0, 0, -t123, -t114 * qJD(3) * t116 - t113 * t122, 0, 0, 0;];
JgD_rot  = t1;
