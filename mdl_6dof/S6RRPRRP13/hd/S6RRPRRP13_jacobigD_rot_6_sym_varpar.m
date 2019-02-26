% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t161 = sin(pkin(6));
t165 = sin(qJ(1));
t181 = t161 * t165;
t167 = cos(qJ(2));
t180 = t161 * t167;
t168 = cos(qJ(1));
t179 = t161 * t168;
t164 = sin(qJ(2));
t178 = t165 * t164;
t177 = t165 * t167;
t176 = t167 * t168;
t175 = t168 * t164;
t174 = qJD(1) * t161;
t166 = cos(qJ(4));
t173 = qJD(2) * t166;
t162 = cos(pkin(6));
t172 = t162 * t176 - t178;
t171 = t162 * t177 + t175;
t170 = t162 * t175 + t177;
t169 = -t162 * t178 + t176;
t163 = sin(qJ(4));
t1 = [0, t168 * t174, 0, -t170 * qJD(1) - t171 * qJD(2) (t171 * t163 + t166 * t181) * qJD(4) - t169 * t173 + (t163 * t179 - t172 * t166) * qJD(1), 0; 0, t165 * t174, 0, t169 * qJD(1) + t172 * qJD(2) (-t172 * t163 - t166 * t179) * qJD(4) - t170 * t173 + (t163 * t181 - t171 * t166) * qJD(1), 0; 0, 0, 0, qJD(2) * t180, -t161 * t164 * t173 + (t162 * t166 - t163 * t180) * qJD(4), 0;];
JgD_rot  = t1;
