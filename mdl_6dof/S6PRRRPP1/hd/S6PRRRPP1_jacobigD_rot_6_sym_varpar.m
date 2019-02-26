% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPP1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:45
% EndTime: 2019-02-26 20:08:45
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t163 = sin(pkin(6));
t166 = sin(qJ(3));
t176 = t163 * t166;
t167 = sin(qJ(2));
t175 = t163 * t167;
t165 = cos(pkin(6));
t174 = t165 * t167;
t169 = cos(qJ(2));
t173 = t165 * t169;
t172 = qJD(2) * t166;
t162 = sin(pkin(10));
t164 = cos(pkin(10));
t171 = t162 * t169 + t164 * t174;
t170 = -t162 * t174 + t164 * t169;
t168 = cos(qJ(3));
t1 = [0, 0, t170 * qJD(2) (t162 * t176 + t170 * t168) * qJD(3) + (-t162 * t173 - t164 * t167) * t172, 0, 0; 0, 0, t171 * qJD(2) (-t164 * t176 + t171 * t168) * qJD(3) + (-t162 * t167 + t164 * t173) * t172, 0, 0; 0, 0, qJD(2) * t175, t163 * t169 * t172 + (t165 * t166 + t168 * t175) * qJD(3), 0, 0;];
JgD_rot  = t1;
