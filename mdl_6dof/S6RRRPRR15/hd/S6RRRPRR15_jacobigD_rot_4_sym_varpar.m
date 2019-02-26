% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR15_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobigD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:20
% EndTime: 2019-02-26 22:24:20
% DurationCPUTime: 0.04s
% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
t166 = sin(pkin(6));
t179 = t166 * cos(pkin(7));
t169 = sin(qJ(2));
t170 = sin(qJ(1));
t178 = t169 * t170;
t172 = cos(qJ(1));
t177 = t169 * t172;
t171 = cos(qJ(2));
t176 = t170 * t171;
t175 = t171 * t172;
t174 = qJD(1) * t166;
t165 = sin(pkin(7));
t173 = qJD(2) * t165;
t168 = cos(pkin(6));
t1 = [0, t172 * t174 -(t168 * t178 - t175) * t173 + (-(-t168 * t175 + t178) * t165 + t172 * t179) * qJD(1), 0, 0, 0; 0, t170 * t174 -(-t168 * t177 - t176) * t173 + (-(-t168 * t176 - t177) * t165 + t170 * t179) * qJD(1), 0, 0, 0; 0, 0, t166 * t169 * t173, 0, 0, 0;];
JgD_rot  = t1;
