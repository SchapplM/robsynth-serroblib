% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRPRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.08s
% Computational Cost: add. (19->12), mult. (73->29), div. (0->0), fcn. (77->12), ass. (0->18)
t168 = sin(pkin(11));
t175 = cos(pkin(6));
t181 = t168 * t175;
t166 = sin(pkin(13));
t171 = cos(pkin(13));
t176 = sin(qJ(3));
t177 = cos(qJ(3));
t178 = qJD(3) * (-t166 * t177 - t171 * t176);
t163 = sin(pkin(7)) * t178;
t170 = sin(pkin(6));
t180 = t170 * t163;
t173 = cos(pkin(11));
t179 = t173 * t175;
t172 = cos(pkin(12));
t167 = sin(pkin(12));
t165 = (t166 * t176 - t171 * t177) * qJD(3);
t164 = cos(pkin(7)) * t178;
t1 = [0, 0, 0, 0 -(-t167 * t181 + t172 * t173) * t165 - (-t167 * t173 - t172 * t181) * t164 - t168 * t180, 0; 0, 0, 0, 0 -(t167 * t179 + t168 * t172) * t165 - (-t167 * t168 + t172 * t179) * t164 + t173 * t180, 0; 0, 0, 0, 0, -t175 * t163 + (-t164 * t172 - t165 * t167) * t170, 0;];
JgD_rot  = t1;
