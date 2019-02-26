% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR9_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:13
% EndTime: 2019-02-26 21:33:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (7->7), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
t106 = sin(qJ(2));
t107 = sin(qJ(1));
t114 = t106 * t107;
t109 = cos(qJ(1));
t113 = t106 * t109;
t108 = cos(qJ(2));
t112 = t107 * t108;
t111 = t108 * t109;
t104 = sin(pkin(6));
t110 = qJD(1) * t104;
t105 = cos(pkin(6));
t1 = [0, t109 * t110, 0, 0 (t105 * t114 - t111) * qJD(2) + (-t105 * t111 + t114) * qJD(1), 0; 0, t107 * t110, 0, 0 (-t105 * t113 - t112) * qJD(2) + (-t105 * t112 - t113) * qJD(1), 0; 0, 0, 0, 0, -t104 * qJD(2) * t106, 0;];
JgD_rot  = t1;
