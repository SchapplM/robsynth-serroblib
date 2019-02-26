% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:21
% EndTime: 2019-02-26 20:16:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t161 = sin(pkin(6));
t164 = sin(qJ(3));
t174 = t161 * t164;
t165 = sin(qJ(2));
t173 = t161 * t165;
t163 = cos(pkin(6));
t172 = t163 * t165;
t167 = cos(qJ(2));
t171 = t163 * t167;
t170 = qJD(2) * t164;
t160 = sin(pkin(11));
t162 = cos(pkin(11));
t169 = t160 * t167 + t162 * t172;
t168 = -t160 * t172 + t162 * t167;
t166 = cos(qJ(3));
t159 = t161 * t167 * t170 + (t163 * t164 + t166 * t173) * qJD(3);
t158 = (t160 * t174 + t168 * t166) * qJD(3) + (-t160 * t171 - t162 * t165) * t170;
t157 = (-t162 * t174 + t169 * t166) * qJD(3) + (-t160 * t165 + t162 * t171) * t170;
t1 = [0, 0, t168 * qJD(2), t158, t158, 0; 0, 0, t169 * qJD(2), t157, t157, 0; 0, 0, qJD(2) * t173, t159, t159, 0;];
JgD_rot  = t1;
