% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR4_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->8), mult. (104->22), div. (0->0), fcn. (112->8), ass. (0->19)
t158 = sin(pkin(12));
t160 = cos(pkin(12));
t162 = sin(qJ(2));
t164 = cos(qJ(2));
t167 = t164 * t158 + t162 * t160;
t169 = qJD(2) * t167;
t159 = sin(pkin(6));
t168 = qJD(1) * t159;
t166 = t158 * t162 - t160 * t164;
t165 = cos(qJ(1));
t163 = sin(qJ(1));
t161 = cos(pkin(6));
t156 = t166 * qJD(2);
t155 = t166 * t161;
t154 = t161 * t169;
t153 = t159 * t169;
t152 = t165 * t154 - t163 * t156 + (-t155 * t163 + t165 * t167) * qJD(1);
t151 = -t163 * t154 - t165 * t156 + (-t155 * t165 - t163 * t167) * qJD(1);
t1 = [0, t165 * t168, 0, t151, t151, 0; 0, t163 * t168, 0, t152, t152, 0; 0, 0, 0, t153, t153, 0;];
JgD_rot  = t1;
