% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:52
% EndTime: 2019-02-26 20:19:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (32->10), mult. (120->29), div. (0->0), fcn. (126->8), ass. (0->19)
t159 = sin(pkin(6));
t162 = sin(qJ(3));
t172 = t159 * t162;
t163 = sin(qJ(2));
t171 = t159 * t163;
t161 = cos(pkin(6));
t170 = t161 * t163;
t165 = cos(qJ(2));
t169 = t161 * t165;
t168 = qJD(2) * t162;
t158 = sin(pkin(12));
t160 = cos(pkin(12));
t167 = t158 * t165 + t160 * t170;
t166 = -t158 * t170 + t160 * t165;
t164 = cos(qJ(3));
t157 = t159 * t165 * t168 + (t161 * t162 + t164 * t171) * qJD(3);
t156 = (t158 * t172 + t166 * t164) * qJD(3) + (-t158 * t169 - t160 * t163) * t168;
t155 = (-t160 * t172 + t167 * t164) * qJD(3) + (-t158 * t163 + t160 * t169) * t168;
t1 = [0, 0, t166 * qJD(2), t156, t156, t156; 0, 0, t167 * qJD(2), t155, t155, t155; 0, 0, qJD(2) * t171, t157, t157, t157;];
JgD_rot  = t1;
