% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:52
% EndTime: 2019-02-26 22:08:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t157 = sin(pkin(6));
t160 = sin(qJ(2));
t177 = t157 * t160;
t161 = sin(qJ(1));
t176 = t157 * t161;
t164 = cos(qJ(1));
t175 = t157 * t164;
t174 = t160 * t161;
t173 = t160 * t164;
t163 = cos(qJ(2));
t172 = t161 * t163;
t171 = t164 * t163;
t170 = qJD(1) * t157;
t162 = cos(qJ(3));
t169 = qJD(2) * t162;
t158 = cos(pkin(6));
t168 = t158 * t171 - t174;
t167 = t158 * t172 + t173;
t166 = t158 * t173 + t172;
t165 = -t158 * t174 + t171;
t159 = sin(qJ(3));
t1 = [0, t164 * t170, t168 * qJD(1) + t165 * qJD(2), 0, 0 (-t165 * t159 + t162 * t176) * qJD(3) - t167 * t169 + (t159 * t175 - t166 * t162) * qJD(1); 0, t161 * t170, t167 * qJD(1) + t166 * qJD(2), 0, 0 (-t166 * t159 - t162 * t175) * qJD(3) + t168 * t169 + (t159 * t176 + t165 * t162) * qJD(1); 0, 0, qJD(2) * t177, 0, 0, t157 * t163 * t169 + (t158 * t162 - t159 * t177) * qJD(3);];
JgD_rot  = t1;
